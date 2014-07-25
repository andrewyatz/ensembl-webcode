=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package Bio::EnsEMBL::ExternalData::BigFile::BigWigAdaptor;
use strict;

use SiteDefs;
use Data::Dumper;
use Bio::DB::BigFile;
use Bio::DB::BigFile::Constants;
use POSIX qw/ceil/;
my $DEBUG = 0;
our $USE_PP = 0;


sub new {
  my ($class, $url) = @_;

  my $self = bless {
    _cache => {},
    _url => $url,
  }, $class;
      
  return $self;
}

sub url { return $_[0]->{'_url'} };

sub bigwig_open {
  my $self = shift;

  Bio::DB::BigFile->set_udc_defaults;
  $self->{_cache}->{_bigwig_handle} ||= Bio::DB::BigFile->bigWigFileOpen($self->url);
  return $self->{_cache}->{_bigwig_handle};
}


# UCSC prepend 'chr' on human chr ids. These are in some of the BigWig
# files. This method returns a possibly modified chr_id after
# checking whats in the BigWig file
sub munge_chr_id {
  my ($self, $chr_id) = @_;
  my $bw = $self->bigwig_open;
  
  warn "Failed to open BigWig file " . $self->url unless $bw;
  
  return undef unless $bw;

  my $list = $bw->chromList;
  my $head = $list->head;
  my $ret_id;
  
  do {
    $ret_id = $head->name if $head->name =~ /^(chr)?$chr_id$/ && $head->size; # Check we get values back for seq region. Maybe need to add 'chr' 
  } while (!$ret_id && ($head = $head->next));
  
  warn " *** could not find region $chr_id in BigWig file" unless $ret_id;
  
  return $ret_id;
}

sub fetch_extended_summary_array {
  my ($self, $chr_id, $start, $end, $bins) = @_;

  my $bw = $self->bigwig_open;
  warn "Failed to open BigWig file" . $self->url unless $bw;
  return [] unless $bw;
  
  #  Maybe need to add 'chr' 
  my $seq_id = $self->munge_chr_id($chr_id);
  return [] if !defined($seq_id);

# Remember this method takes half-open coords (subtract 1 from start)
  my $summary_e = $bw->bigWigSummaryArrayExtended("$seq_id",$start-1,$end,$bins);

  if ($DEBUG) {
    warn " *** fetch extended summary: $chr_id:$start-$end : found ", scalar(@$summary_e), " summary points\n";
  }
  
  return $summary_e;
}

=head2 _pp_mean_summary
 
  Args [1]   : Bio::EnsEMBL::Slice; Slice to query against
  Args [2]   : Integer; number of summary bins required
  Args [3]   : Bio::EnsEMBL::Analysis; Analysis object to attach to the eventual data hash. Can be undef
  Description: Retrieve a summary of the requested region as mean scores within the region
  Returntype : ArrayRef of summary hashes. Keys are min, max, features. Features 
               is an ArrayRef holding the keys start, end, slice, analysis, score representing
               the summarised averaged block
  Exceptions : None

=cut


sub mean_summary {
  my ($self, $slice, $bins, $analysis) = @_;
  my $bw = $self->bigwig_open;
  warn "Failed to open BigWig file" . $self->url unless $bw;
  return { features => [], min => 0, max => 0 } unless $bw;
  
  #  Maybe need to add 'chr' 
  my $seq_id = $self->munge_chr_id($slice->seq_region_name());
  return { features => [], min => 0, max => 0 } if !defined($seq_id);

  # Query
  my $summary;
  my @params = ($bw, $seq_id, $slice->start-1, $slice->end, $bins, $slice->strand, $slice, $analysis);
  if($USE_PP) { # if forced use Perl
    $summary = $self->_pp_mean_summary(@params);
  }
  else { # if not use what's available
    $summary = $self->_mean_summary(@params);
  }
  if ($DEBUG) {
    warn " *** fetch extended summary: ".$slice->name()." : found ", scalar(@{$summary->{features}}), " summary points\n";
  }

  return $summary;
}

=head2 _pp_mean_summary
 
  Args [1]   : Bio::DB::BigFile; BigWig file to query against
  Args [2]   : String; seq_region to query for
  Args [3]   : Integer; start to query (UCSC formatted)
  Args [4]   : Integer; end to query to (UCSC formatted)
  Args [5]   : Integer; number of summary bins required
  Args [6]   : Integer; strand to produce elements for (1/-1)
  Args [7]   : Bio::EnsEMBL::Slice; Slice object to attach to the eventual data hash. Can be undef
  Args [8]   : Bio::EnsEMBL::Analysis; Analysis object to attach to the eventual data hash. Can be undef
  Description: A pure perl implementation to summarise a BigWig file into mean bins of data. This
               method is used only when Inline::C is not available to compile a working C based
               solution (which is faster than this).
  Returntype : ArrayRef of summary hashes. Keys are min, max, features. Features 
               is an ArrayRef holding the keys start, end, slice, analysis, score representing
               the summarised averaged block
  Exceptions : None

=cut

sub _pp_mean_summary {
  my ($self, $bw_file, $seq_region, $start, $end, $bins, $strand, $slice, $analysis) = @_;
  my $summary = $bw_file->bigWigSummaryArray($seq_region, $start, $end, bbiSumMean, $bins);
  my $query_width = ($end-$start);
  my $bin_width = ceil(($query_width) / $bins);
  my $flip      = $strand == -1 ? $query_width + 1 : undef;
  my @features;
  
  my $min_max_init = 0;
  my ($min, $max) = (0,0);
  for (my $i = 0; $i < $bins; $i++) {
    my $score = $summary->[$i];
    next unless defined $score;

    if(!$min_max_init) {
      ($min, $max) = ($score, $score);
      $min_max_init = 1;
    }

    $min = $score if $score < $min;
    $max = $score if $score > $max;
    
    push @features, {
      start => $flip ? $flip - (($i + 1) * $bin_width) : ($i * $bin_width + 1),
      end   => $flip ? $flip - ($i * $bin_width + 1)   : (($i + 1) * $bin_width),
      score => $score,
      analysis => $analysis,
      slice => $slice
    };
  }

  return {
    min => $min,
    max => $max,
    features => \@features
  };
}

eval {
   require Inline;
   Inline->import (C => Config =>
                   BUILD_NOISY => 1,
                   INC => "-I$SiteDefs::JKSRC_DIR/inc",
                   LIBS => "-L$SiteDefs::JKSRC_DIR/lib/$ENV{MACHTYPE}/jkweb.a",
                   DIRECTORY => "$SiteDefs::ENSEMBL_WEBROOT/cbuild",
                   TYPEMAPS => "$SiteDefs::ENSEMBL_WEBROOT/typemap");
   Inline->import (C =><<'EOC');

#include "common.h"
#include "bbiFile.h"
#include "bigWig.h"
#include "math.h"

typedef struct bbiFile     *Bio__DB__bbiFile;

SV * _mean_summary(SV* self, SV* bwf, char* seq_region, int start, int end, int bins, int strand, SV* slice, SV* analysis) {
  int     i;
  int     binWidth;
  int     queryWidth;
  int     summaryType;
  double  min;
  double  max;
  double  currentValue;
  boolean result;
  double  *values;
  HV     *h;
  AV     *av;

  queryWidth = (end - start);
  binWidth = ceil((double)queryWidth/(double)bins);
  summaryType = 0;

  //We need to access the backing C struct not the Perl SV* ref
  struct bbiFile *unwrappedBbiFile;
  unwrappedBbiFile = (struct bbiFile*)SvIV(SvRV(bwf));

  //values is the elements brought back from BigWig according to the summary
  values = Newx(values,bins,double);
  //This is the call to kent src libs to get the summary into values
  result = bigWigSummaryArray(unwrappedBbiFile,seq_region,start,end,summaryType,bins,values);

  // Get our return values ready
  av = newAV();
  min = 0.0;
  max = 0.0;
   
  if (result == TRUE) {
   
   for (i=0;i<bins;i++) {
     currentValue = values[i];
     if(! isnan(currentValue)) {
       h = newHV();
       //Check for the current min value
       if(min == 0) {
         min = currentValue;
       }
       else {
         if(currentValue < min) {
           min = currentValue;
         }
       }
       
       //Do the same for max
       if(max == 0) {
         max = currentValue;
       }
       else {
         if(currentValue > max) {
           max = currentValue;
         }
       }
       
       // Calculate bin coordinates in Ensembl positions (+1)
       unsigned int start,end;
       if(strand == 1) {
         start = i*binWidth+1; //ensembl start
         end = (i+1)*binWidth; //ensembl end
       }
       else {
         start = queryWidth - ((i+1)*binWidth); //ensembl start
         end = queryWidth - (i*binWidth+1); //ensembl end
       }
       
       //Now store them in the hash and add the hash to the output array
       (void)hv_store(h,"score",    5, newSVnv(currentValue),0);
       (void)hv_store(h,"start",    5, newSVnv(start),0);
       (void)hv_store(h,"end",      3, newSVnv(end),0);
       (void)hv_store(h,"slice",    5, (!SvOK(slice) ? &PL_sv_undef : SvREFCNT_inc(slice)),0);
       (void)hv_store(h,"analysis", 8, (!SvOK(analysis) ? &PL_sv_undef : SvREFCNT_inc(analysis)),0);
       av_push(av, newRV_noinc((SV*)h));
     }
   }
 }
 Safefree(values);
 HV *returnHash = newHV();
 (void)hv_store(returnHash, "min", 3, newSVnv(min),0);
 (void)hv_store(returnHash, "max", 3, newSVnv(max),0);
 (void)hv_store(returnHash, "features", 8, (SV*) newRV_noinc((SV*)av),0);
 return (SV*) newRV_noinc((SV*)returnHash);
}

EOC
};

if($@) {
  warn $@;
  *_mean_summary =\&_pp_mean_summary;
}


1;
