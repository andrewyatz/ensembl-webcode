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

package EnsEMBL::Draw::GlyphSet::regulatory_regions;

use strict;

use base qw(EnsEMBL::Draw::GlyphSet_simple);

sub squish { return 1; }

sub get_feature_sets {
  my ($self, $fg_db) = @_;

  my $logic_name = $self->my_config('logic_name')
                   || $self->my_config('description');
  my $fg_a_a =  $fg_db->get_AnalysisAdaptor;
  my $fg_fs_a = $fg_db->get_FeatureSetAdaptor;
  my $analysis = $fg_a_a->fetch_by_logic_name($logic_name);
  return [grep {
    $_->name !~ /cisRED\s+search\s+regions/i
  } @{$fg_fs_a->fetch_all_by_feature_class('external',undef,{
    constraints => {
      analyses => [$analysis],
    },
  })}];
}

sub features {
  my ($self) = @_;
  my $slice = $self->{'container'};
  my $wuc = $self->{'config'};
 
  my $efg_db = undef;
  my $db_type  = $self->my_config('db_type')||'funcgen';
  unless($slice->isa("Bio::EnsEMBL::Compara::AlignSlice::Slice")) {
    $efg_db = $slice->adaptor->db->get_db_adaptor($db_type);
    if(!$efg_db) {
      warn("Cannot connect to $db_type db");
      return [];
    }
  }


  my @fsets = @{$self->get_feature_sets($efg_db)}; 
  my $external_Feature_adaptor  = $efg_db->get_ExternalFeatureAdaptor;
  my $f = $external_Feature_adaptor->fetch_all_by_Slice_FeatureSets($slice, \@fsets);

  # count used for colour assignment
  my $count = 0;
  foreach my $feat (@$f){
    $wuc->cache($feat->display_label, $count);   
    $count ++;
    if ($count >= 15) {$count = 0;} 
  } 

  return $f;
}

sub href {
  my ($self, $f) = @_;
  my $id = $f->display_label;
  my $dbid = $f->dbID;
  my $analysis =  $f->analysis->logic_name;

  my $href = $self->_url
  ({'action'   => 'RegFeature',
    'fid'      => $id,
    'ftype'    => $analysis,
    'dbid'     => $dbid, 
    'species'  => $self->species, 
  });

  return $href;
}



sub colour_key {
  my ($self, $f) = @_;
  my $wuc = $self->{'config'}; 
  my $colour = $wuc->cache($f->display_label); 
  return $colour;
}



1;
