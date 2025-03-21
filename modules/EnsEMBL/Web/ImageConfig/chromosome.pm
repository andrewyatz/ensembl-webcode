=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

package EnsEMBL::Web::ImageConfig::chromosome;

use strict;
use warnings;

use parent qw(EnsEMBL::Web::ImageConfig);

sub init_extra_menus {
  shift->add_extra_menu('display_option');
}

sub init_cacheable {
  ## @override
  my $self = shift;

  $self->SUPER::init_cacheable(@_);

  $self->create_menus('decorations');

  $self->add_tracks('decorations',
    [ 'ideogram', 'Ideogram', 'ideogram',  { display => 'normal', menu => 'no', strand => 'r', colourset => 'ideogram' }],
  );

  $self->load_tracks;

  $self->add_tracks('decorations',
    [ 'draggable', '', 'draggable', { display => 'normal', menu => 'no' }]
  );

  $self->get_node('decorations')->set_data('caption', 'Decorations');

  $self->modify_configs(
    [ 'decorations' ],
    { short_labels => 1 }
  );
}

1;
