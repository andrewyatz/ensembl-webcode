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

package Bio::EnsEMBL::GlyphSet::P_domain;

use strict;

use base qw(Bio::EnsEMBL::GlyphSet);

## Variables defined in UserConfig.pm 
## 'caption'   -> Track label
## 'logicname' -> Logic name

sub colour_key { return $_[1]->analysis->logic_name; }

sub _init {
  my $self = shift;
  return $self->render_text if $self->{'text_export'};
  my $protein = $self->{'container'};

  $self->_init_bump;

  my $label        = $self->my_config('caption');
  my $h            = $self->my_config('height') || 4;
  my $depth        = $self->depth;
  my $font_details = $self->get_text_simple(undef, 'innertext');
  my $pix_per_bp   = $self->scalex;
  my $count = 0;

  foreach my $logic_name (@{$self->my_config('logic_names') || []}) {
    my (%hash, $colour);
    my $prot_features = $protein->get_all_ProteinFeatures($logic_name);
    $count += scalar(@{$prot_features||[]});
    
    push @{$hash{$_->hseqname}}, $_ for @{$prot_features||[]};
    
    foreach my $key (keys %hash) {
      my (@rect, $prsave, $minx, $maxx);
      
      foreach my $pr (@{$hash{$key}}) {
        my $x    = $pr->start;
           $minx = $x if $x < $minx || !defined $minx;
        my $w    = $pr->end - $x;
           $maxx = $pr->end if $pr->end > $maxx || !defined $maxx;
        my $id   = $pr->hseqname;
        
        push @rect, $self->Rect({
          x      => $x,
          y      => 0,
          width  => $w,
          height => $h,
          colour => $colour ||= $self->get_colour($pr),
        });
        
        $prsave ||= $pr;
      }
      
      my $title  = "$label $key; Positions: $minx-$maxx";
         $title .= '; Interpro: ' . $prsave->interpro_ac if $prsave->interpro_ac;
         $title .= '; '. $prsave->idesc                  if $prsave->idesc;
      my $dbID   = $prsave->dbID;
      
      my $composite = $self->Composite({
        x     => $minx,
        y     => 0,
        href  => $self->_url({ type => 'Transcript', action => 'ProteinSummary', pf_id => $dbID }),
        title => $title
      });
      
      $composite->push(@rect,
        $self->Rect({
          x         => $minx,
          y         => $h / 2,
          width     => $maxx - $minx,
          height    => 0,
          colour    => $colour,
          absolutey => 1,
        })
      );
      
      #### add a label
      my $desc = $prsave->idesc || $key;
      my @res  = $self->get_text_width(0, $desc, '', font => $font_details->{'font'}, ptsize => $font_details->{'fontsize'});
      
      $composite->push($self->Text({
        font      => $font_details->{'font'},
        ptsize    => $font_details->{'fontsize'},
        halign    => 'left',
        text      => $desc,
        x         => $composite->x,
        y         => $h,
        height    => $font_details->{'height'},
        width     => $res[2] / $pix_per_bp,
        colour    => $colour,
        absolutey => 1
      }));
      
      if ($depth > 0) {
        my $bump_start = int($composite->x * $pix_per_bp);
        my $bump_end   = $bump_start + int($composite->width / $pix_per_bp);
        my $row        = $self->bump_row($bump_start, $bump_end);
        
        $composite->y($composite->y + ($row * (4 + $h + $font_details->{'height'}))) if $row;
      }
      
      $self->push($composite);
    }
  }
  $self->no_features unless $count;
}

sub render_text {
  my $self      = shift;
  my $container = $self->{'container'};
  my $label     = $self->my_config('caption');
  my $export;
  
  foreach my $logic_name (@{$self->my_config('logic_names') || []}) {
    my @features = map { $_->[1] } sort { $a->[0] cmp $b->[0] } map { [ $_->hseqname, $_ ] } @{$container->get_all_ProteinFeatures($logic_name)};
    
    foreach (@features) {
      my $analysis = $_->analysis;
      
      $export .= $self->_render_text($_, $analysis->gff_feature, { 
        headers => [ 'id', 'description' ],
        values  => [ $_->hseqname, $_->idesc ]
      }, {
        source  => $analysis->gff_source,
      });
    }
  }
  
  return $export;
}

1;
