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


package EnsEMBL::Web::Component::StructuralVariation::SupportingEvidence;

use strict;

use base qw(EnsEMBL::Web::Component);

sub _init {
  my $self = shift;
  $self->cacheable(0);
  $self->ajaxable(1);
}

sub content {
  my $self = shift;
  my $object        = $self->object;
  my $hub           = $self->hub;
  my $supporting_sv = $object->supporting_sv;
  my $html          = $self->supporting_evidence_table($supporting_sv);
  return $html;
}


sub supporting_evidence_table {
  my $self     = shift;
  my $ssvs     = shift;
  my $hub      = $self->hub;
  my $object   = $self->object;
  my $title    = 'Supporting evidence';
  my $table_id = 'evidence';
  
  my $columns = [
     { key => 'ssv',    sort => 'string',        title => 'Supporting evidence'   },
     { key => 'pos',    sort => 'position_html', title => 'Chr:bp (strand)'       }, 
  ];
  my $sorting = 'pos';

  my $rows = ();
  my $has_bp;
  my $has_phen;
  
  # Supporting evidences list
  if (scalar @{$ssvs}) {
    my $ssv_names = {};
    foreach my $ssv (@$ssvs){
      my $name = $ssv->variation_name;
      $name =~ /(\d+)$/;
      my $ssv_nb = $1;
      $ssv_names->{$1} = $ssv;
    }
   
    foreach my $ssv_n (sort {$a <=> $b} (keys(%$ssv_names))) {
      my $ssv_obj = $ssv_names->{$ssv_n};
      my $name = $ssv_obj->variation_name;
      my $loc;
      my $bp_order;
      my $is_somatic = $ssv_obj->is_somatic;
      
      # Location(s)
      foreach my $svf (sort {$a->seq_region_start <=> $b->seq_region_start} @{$ssv_obj->get_all_StructuralVariationFeatures}) {
        my ($sv_name,$chr_bp);
        my $start  = $svf->seq_region_start;
        my $end    = $svf->seq_region_end;
        my $strand = $svf->seq_region_strand;
           $strand = ' <small>('.($strand > 0 ? 'forward' : 'reverse').' strand)</small>';
        next if ($start == 0 || $end == 0);
        
        $bp_order = $svf->breakpoint_order;
        my $chr = $svf->seq_region_name.':';
        
        
        if ($bp_order && $is_somatic == 1) {
          my @c_list = ($start!=$end) ? ($start,$end) : ($start);
            
          foreach my $coord (@c_list) {
            $chr_bp = $chr.$coord;
            my $loc_hash = {
                type   => 'Location',
                action => 'View',
                r      => $chr_bp,
            };

            $loc_hash->{sv} .= $name if ($ssv_obj->is_evidence == 0);
          
            my $loc_url = $hub->url($loc_hash);
            $loc .= ($loc) ? "<br />to " : "from ";
            $loc .= qq{<a href="$loc_url">$chr_bp</a>$strand};
          }
        }
        else {
          $chr_bp  = $chr;
          $chr_bp .= $start == $end ? $start : $start.'-'.$end;
          my $loc_url;
        
          my $loc_hash = {
              type   => 'Location',
              action => 'View',
              r      => $chr_bp,
          };
          $loc_hash->{sv} = $name if ($ssv_obj->is_evidence == 0);
          
          my $loc_url = $hub->url($loc_hash);
          $loc .= "<br />" if ($loc);
          $loc .= qq{<a href="$loc_url">$chr_bp</a>$strand};
        }
        
        $has_bp = 1 if ($bp_order && $bp_order > 1 && $is_somatic == 1);
      }
      $loc = '-' if (!$loc);
  
      # Name
      if ($ssv_obj->is_evidence == 0) {
        my $sv_link = $hub->url({
                        type   => 'StructuralVariation',
                        action => 'Explore',
                        sv     => $name,
                      });
        $name = qq{<a href="$sv_link">$name</a>};
      }
  
      # Class + class colour
      my $colour = $object->get_class_colour($ssv_obj->class_SO_term);
      my $sv_class = '<div><div style="float:left;background-color:'.$colour.';padding:5px;margin-top:4px"></div> <div style="float:left;margin-left:5px">'.$ssv_obj->var_class.'</div></div>';
       
      # Annotation(s)
      my $clin = $ssv_obj->clinical_significance;
      my ($indiv, $strain, $phen);
      my ($indivs, $strains, $phens);
      
      foreach my $pf (sort {$a->seq_region_start <=> $b->seq_region_start} @{$ssv_obj->get_all_PhenotypeFeatures}) {
        my $a_phen = $pf->phenotype->description;
        $phens->{$a_phen} = 1;
        $has_phen = 1;
      }
      
      foreach my $svs (@{$ssv_obj->get_all_StructuralVariationSamples}) {
        
        my $a_indiv  = ($svs->individual) ? $svs->individual->name : undef;
        my $a_strain = ($svs->strain) ? $svs->strain->name : undef;
        
        $indivs->{$a_indiv} = 1 if ($a_indiv);
        $strains->{$a_strain} = 1 if ($a_strain); 
      }
      
      $indiv  = join(';<br />', sort (keys(%$indivs)));
      $strain = join(';<br />', sort (keys(%$strains)));
      $phen   = join(';<br />', sort (keys(%$phens)));
      
      my %row = (
                  ssv      => $name,
                  class    => $sv_class,
                  pos      => $loc,
                  clin     => $clin ? $clin : '-',
                  indiv    => $indiv ? $indiv : '-',
                  strain   => $strain ? $strain : '-',
                  bp_order => $bp_order ? $bp_order : '-',
                  phen     => $phen ? $phen : '-',
                );
        
      push @$rows, \%row;
    }
    
    if (defined($has_bp)) {  
     push(@$columns,{ key => 'bp_order', sort => 'numeric', title => 'Breakpoint order' });
     $sorting = 'bp_order';
    };
    
    push @$columns, (
     { key => 'class', sort => 'string', title => 'Allele type'           },
     { key => 'clin',  sort => 'string', title => 'Clinical significance' },
     { key => 'indiv', sort => 'string', title => 'Individual name(s)'    },
    );
    
    if (defined($has_phen)) {  
     push(@$columns,{ key => 'phen', sort => 'string', title => 'Phenotype(s)' });
    };
    
    if ($self->hub->species ne 'Homo_sapiens') {  
     push(@$columns,{ key => 'strain', sort => 'string', title => 'Strain' });
    };
    
    
    return $self->new_table($columns, $rows, { data_table => 1, sorting => [ "$sorting asc" ] })->render;
  }
}
1;
