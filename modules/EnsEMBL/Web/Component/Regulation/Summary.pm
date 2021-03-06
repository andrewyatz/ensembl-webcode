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

package EnsEMBL::Web::Component::Regulation::Summary;

use strict;
use warnings;
no warnings "uninitialized";
use base qw(EnsEMBL::Web::Component::Regulation);


sub _init {
  my $self = shift;
  $self->cacheable( 0 );
  $self->ajaxable(  0 );
}

sub content {
  my $self    = shift;
  my $object  = $self->object;
  my $summary = $self->new_twocol;

  my $url = $self->hub->url({
    'type'   => 'Location',
    'action' => 'View',
    'r'      => $object->seq_region_name.':'.$object->bound_start.'-'.$object->bound_end
  });

  my $location_html = sprintf('<p><a href="%s" class="constant">%s: %s-%s</a></p>',
    $url,
    $object->neat_sr_name( $object->seq_region_type, $object->seq_region_name ),
    $object->thousandify( $object->seq_region_start ),
    $object->thousandify( $object->seq_region_end ),
  );

  $summary->add_row('Location', $location_html);

  ## TODO - use a simpler method to create simple tables
  my $table = $self->dom->create_element('table', {'cellspacing' => '0', 'children' => [
    {'node_name' => 'tr', 'children' => [ map {'node_name' => 'th', 'inner_text' => $_}, ('Cell line', 'Feature type', 'Bound co-ordinates') ]}
  ]});

  my $all_objs = $object->fetch_all_objs;
  foreach my $reg_object (sort { $a->feature_set->cell_type->name cmp $b->feature_set->cell_type->name } @$all_objs ) {
    next if $reg_object->feature_set->cell_type->name =~/MultiCell/;
    $table->append_child('tr', {'children' => [map {'node_name' => 'td', 'style' => 'padding-right: 8px', 'inner_text' => $_},
      $reg_object->feature_set->cell_type->name,
      $reg_object->feature_type->name,
      $reg_object->bound_start ."-". $reg_object->bound_end
    ]});
  }

  $summary->add_row('', $table->render);

  return $summary->render;
}

1;
