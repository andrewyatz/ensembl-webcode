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

package EnsEMBL::Web::ImageConfig;

use strict;

use Digest::MD5 qw(md5_hex);
use HTML::Entities qw(encode_entities decode_entities);
use JSON qw(from_json);
use URI::Escape qw(uri_unescape);

use Bio::EnsEMBL::ExternalData::DAS::Coordinator;

use Sanger::Graphics::TextHelper;

use EnsEMBL::Web::DBSQL::DBConnection;
use EnsEMBL::Web::Tree;

#########
# 'user' settings are restored from cookie if available
#  default settings are overridden by 'user' settings
#

# Takes two parameters
# (1) - the hub (i.e. an EnsEMBL::Web::Hub object)
# (2) - the species to use (defaults to the current species)

sub new {
  my $class   = shift;
  my $hub     = shift;
  my $species = shift || $hub->species;
  my $code    = shift;
  my $type    = $class =~ /([^:]+)$/ ? $1 : $class;
  my $style   = $hub->species_defs->ENSEMBL_STYLE || {};
  
  my $self = {
    hub              => $hub,
    _font_face       => $style->{'GRAPHIC_FONT'} || 'Arial',
    _font_size       => ($style->{'GRAPHIC_FONTSIZE'} * $style->{'GRAPHIC_LABEL'}) || 20,
    _texthelper      => Sanger::Graphics::TextHelper->new,
    code             => $code,
    type             => $type,
    species          => $species,
    altered          => 0,
    _tree            => EnsEMBL::Web::Tree->new,
    transcript_types => [qw(transcript alignslice_transcript tsv_transcript gsv_transcript TSE_transcript)],
    _parameters      => { # Default parameters
      storable     => 1,      
      has_das      => 1,
      datahubs     => 0,
      image_width  => $ENV{'ENSEMBL_IMAGE_WIDTH'} || 800,
      image_resize => 0,      
      margin       => 5,
      spacing      => 2,
      label_width  => 113,
      show_labels  => 'yes',
      slice_number => '1|1',
      toolbars     => { top => 1, bottom => 0 },
    },
    extra_menus => {
      active_tracks    => 1,
      favourite_tracks => 1,
      search_results   => 1,
      display_options  => 1,
    },
    unsortable_menus => {
      decorations => 1,
      information => 1,
      options     => 1,
      other       => 1,
    },
    alignment_renderers => [
      'off',         'Off',
      'normal',      'Normal',
      'labels',      'Labels',
      'half_height', 'Half height',
      'stack',       'Stacked',
      'unlimited',   'Stacked unlimited',
      'ungrouped',   'Ungrouped',
    ],
  };
  
  return bless $self, $class;
}

sub initialize {
  my $self      = shift;
  my $class     = ref $self;
  my $species   = $self->species;
  my $code      = $self->code;
  my $cache     = $self->hub->cache;
  my $cache_key = $self->cache_key;
  
  # Check memcached for defaults
  if (my $defaults = $cache && $cache_key ? $cache->get($cache_key) : undef) {
    my $user_data = $self->tree->user_data;
    
    $self->{$_} = $defaults->{$_} for keys %$defaults;
    $self->tree->push_user_data_through_tree($user_data);
  } else {
    # No cached defaults found, so initialize them and cache
    $self->init;
    $self->modify;
    
    if ($cache && $cache_key) {
      $self->tree->hide_user_data;
      
      my $defaults = {
        _tree       => $self->{'_tree'},
        _parameters => $self->{'_parameters'},
        extra_menus => $self->{'extra_menus'},
      };
      
      $cache->set($cache_key, $defaults, undef, 'IMAGE_CONFIG', $species);
      $self->tree->reveal_user_data;
    }
  }
  
  my $sortable = $self->get_parameter('sortable_tracks');
  
  $self->set_parameter('sortable_tracks', 1) if $sortable eq 'drag' && $ENV{'HTTP_USER_AGENT'} =~ /MSIE (\d+)/ && $1 < 7; # No sortable tracks on images for IE6 and lower
  $self->{'extra_menus'}{'track_order'} = 1  if $sortable;
  
  $self->{'no_image_frame'} = 1;
  
  # Add user defined data sources
  $self->load_user_tracks;
  
  # Combine info and decorations into a single menu
  my $decorations = $self->get_node('decorations') || $self->get_node('other');
  my $information = $self->get_node('information');
  
  if ($decorations && $information) {
    $decorations->set('caption', 'Information and decorations');
    $decorations->append_children($information->nodes);
  }
}

sub menus {
  return $_[0]->{'menus'} ||= {
    # Sequence
    seq_assembly        => 'Sequence and assembly',
    sequence            => [ 'Sequence',                'seq_assembly' ],
    misc_feature        => [ 'Clones & misc. regions',  'seq_assembly' ],
    genome_attribs      => [ 'Genome attributes',       'seq_assembly' ],
    marker              => [ 'Markers',                 'seq_assembly' ],
    simple              => [ 'Simple features',         'seq_assembly' ],
    ditag               => [ 'Ditag features',          'seq_assembly' ],
    dna_align_other     => [ 'GRC alignments',          'seq_assembly' ],
    dna_align_compara   => [ 'Imported alignments',     'seq_assembly' ],
    
    # Transcripts/Genes
    gene_transcript     => 'Genes and transcripts',
    transcript          => [ 'Genes',                  'gene_transcript' ],
    prediction          => [ 'Prediction transcripts', 'gene_transcript' ],
    Lrg                 => [ 'LRG transcripts',        'gene_transcript' ],
    rnaseq              => [ 'RNASeq models',          'gene_transcript' ],
    
    # Supporting evidence
    splice_sites        => 'Splice sites',
    evidence            => 'Evidence',
    
    # Alignments
    mrna_prot           => 'mRNA and protein alignments',
    dna_align_cdna      => [ 'mRNA alignments',    'mrna_prot' ],
    dna_align_est       => [ 'EST alignments',     'mrna_prot' ],
    protein_align       => [ 'Protein alignments', 'mrna_prot' ],
    protein_feature     => [ 'Protein features',   'mrna_prot' ],
    dna_align_rna       => 'ncRNA',
    
    # Proteins
    domain              => 'Protein domains',
    gsv_domain          => 'Protein domains',
    feature             => 'Protein features',
    
    # Variations
    variation           => 'Variation',
    recombination       => [ 'Recombination & Accessibility', 'variation' ],
    somatic             => 'Somatic mutations',
    ld_population       => 'Population features',
    
    # Regulation
    functional          => 'Regulation',
    
    # Compara
    compara             => 'Comparative genomics',
    pairwise_blastz     => [ 'BLASTz/LASTz alignments',    'compara' ],
    pairwise_other      => [ 'Pairwise alignment',         'compara' ],
    pairwise_tblat      => [ 'Translated blat alignments', 'compara' ],
    multiple_align      => [ 'Multiple alignments',        'compara' ],
    conservation        => [ 'Conservation regions',       'compara' ],
    synteny             => 'Synteny',
    
    # Other features
    repeat              => 'Repeat regions',
    oligo               => 'Oligo probes',
    trans_associated    => 'Transcript features',
    
    # Info/decorations
    information         => 'Information',
    decorations         => 'Additional decorations',
    other               => 'Additional decorations',
    
    # External data
    user_data           => 'Your data',
    external_data       => 'External data',
  };
}

sub init   {}
sub modify {} # For plugins

sub storable     :lvalue { $_[0]{'_parameters'}{'storable'};     } # Set to 1 if configuration can be altered
sub image_resize :lvalue { $_[0]{'_parameters'}{'image_resize'}; } # Set to 1 if there is image resize function
sub has_das      :lvalue { $_[0]{'_parameters'}{'has_das'};      } # Set to 1 if there are DAS tracks
sub altered      :lvalue { $_[0]{'altered'};                     } # Set to 1 if the configuration has been updated

sub hub                 { return $_[0]->{'hub'};                                               }
sub code                { return $_[0]->{'code'};                                              }
sub core_objects        { return $_[0]->hub->core_objects;                                     }
sub colourmap           { return $_[0]->hub->colourmap;                                        }
sub species_defs        { return $_[0]->hub->species_defs;                                     }
sub sd_call             { return $_[0]->species_defs->get_config($_[0]->{'species'}, $_[1]);   }
sub databases           { return $_[0]->sd_call('databases');                                  }
sub texthelper          { return $_[0]->{'_texthelper'};                                       }
sub transform           { return $_[0]->{'transform'};                                         }
sub tree                { return $_[0]->{'_tree'};                                             }
sub species             { return $_[0]->{'species'};                                           }
sub multi_species       { return 0;                                                            }
sub cache_key           { return join '::', '', ref $_[0], $_[0]->species, $_[0]->code;        }
sub bgcolor             { return $_[0]->get_parameter('bgcolor') || 'background1';             }
sub bgcolour            { return $_[0]->bgcolor;                                               }
sub get_node            { return shift->tree->get_node(@_);                                    }
sub get_parameters      { return $_[0]->{'_parameters'};                                       }
sub get_parameter       { return $_[0]->{'_parameters'}{$_[1]};                                }
sub set_width           { $_[0]->set_parameter('width', $_[1]);                                } # TODO: find out why we have width and image_width. delete?
sub image_height        { return shift->parameter('image_height',    @_);                      }
sub image_width         { return shift->parameter('image_width',     @_);                      }
sub container_width     { return shift->parameter('container_width', @_);                      }
sub toolbars            { return shift->parameter('toolbars',        @_);                      }
sub slice_number        { return shift->parameter('slice_number',    @_);                      } # TODO: delete?
sub get_tracks          { return grep { $_->{'data'}{'node_type'} eq 'track' } $_[0]->tree->nodes; } # return a list of track nodes
sub get_sortable_tracks { return grep { $_->get('sortable') && $_->get('menu') ne 'no' } @{$_[0]->glyphset_configs}; }

sub get_user_settings {
  my $self     = shift;
  my $settings = $self->tree->user_data;
  delete $settings->{'track_order'} if $settings->{'track_order'} && !scalar keys %{$settings->{'track_order'}};
  return $settings;
}

sub default_track_order {
  my ($self, @tracks) = @_;
  my (@out, %seen, %before, %after);
  
  foreach (@tracks) {
    next unless my $ta = $_->get('track_after');
    $before{$ta->id} = $_;
    $after{$_->id}   = $ta;
  }
  
  foreach my $t (grep !$after{$_->id}, @tracks) {
    for (my $t2 = $t; $t2 && !$seen{$t2}; $t2 = $before{$t2->id}) {
      push @out, $t2;
      $seen{$t2} = 1; # %seen is paranoia against infinite loops
    }
  }
  
  return @out;
}

sub glyphset_configs {
  my $self = shift;
  
  if (!$self->{'ordered_tracks'}) {
    my @tracks      = $self->get_tracks;
    my $track_order = $self->track_order;
    my $i           = 1;
    my @default_order;
    
    foreach my $track ($self->default_track_order(@tracks)) {
      my $strand = $track->get('strand');
      
      if ($strand =~ /^[rf]$/) {
        if ($strand eq 'f') {
          unshift @default_order, $track;
        } else {
          push @default_order, $track;
        }
      } else {
        my $clone = $self->tree->create_node($track->id . '_tmp'); # Id is replaced in the for loop below
        $clone->{$_} = $_ eq 'data' ? { %{$track->{'data'}} } : $track->{$_} for keys %$track; # Make a new hash for data, so that drawing_strand can differ
        
        $clone->set('drawing_strand', 'f');
        $track->set('drawing_strand', 'r');
        
        unshift @default_order, $clone;
        push    @default_order, $track;
      }
    }
    
    if ($self->get_parameter('sortable_tracks')) {
      $_->set('sortable', 1) for grep !$self->{'unsortable_menus'}->{$_->parent_key}, @default_order;
    }
    
    my @ordered_tracks;
    my %order     = map { join('.', grep $_, $_->id, $_->get('drawing_strand')) => [ $i++, $_ ] } @default_order;
    $order{$_}[0] = $track_order->{$_} for grep $order{$_}, keys %$track_order;
    
    foreach (sort { $a->[0] <=> $b->[0] } values %order) {
      $_->[1]->set('order', $_->[0]);
      push @ordered_tracks, $_->[1];
    }
    
    $self->{'ordered_tracks'} = \@ordered_tracks;
  }
  
  return $self->{'ordered_tracks'};
}

sub get_favourite_tracks {
  my $self = shift;
  
  if (!$self->{'favourite_tracks'}) {
    my $hub  = $self->hub;
    my $user = $hub->user;
    my $data = $hub->session->get_data(type => 'favourite_tracks', code => 'favourite_tracks') || {};
    
    $self->{'favourite_tracks'} = $data->{'tracks'} || {};
    $self->{'favourite_tracks'} = { %{$self->{'favourite_tracks'}}, %{$user->get_favourite_tracks} } if $user;
  }
  
  return $self->{'favourite_tracks'};
}

sub track_order {
  my $self = shift;
  my $node = $self->get_node('track_order');
  return $node ? $node->get($self->species) || {} : {};
}

sub set_user_settings {
  my ($self, $data) = @_;
  
  foreach my $key (keys %$data) {
    my $node = $self->get_node($key) || $self->tree->create_node($key); # If node doesn't exist, it is a track from a different species whose configuration needs to be retained
    
    next unless $node;
    
    my $renderers = $node->data->{'renderers'} || [];
    my %valid     = @$renderers;
    
    foreach (keys %{$data->{$key}}) {
      if ($_ eq 'display' && %valid && !$valid{$data->{$key}{$_}}) {
        $node->set_user($_, $valid{'normal'} ? 'normal' : $renderers->[2]); # index 2 contains the code for the first "on" renderer
      } else {
        $node->set_user($_, $data->{$key}{$_});
      }
    }
  }
}

sub set_parameters {
  my ($self, $params) = @_;
  $self->{'_parameters'}{$_} = $params->{$_} for keys %$params; 
}

sub set_parameter {
  my ($self, $key, $value) = @_;
  $self->{'_parameters'}{$key} = $value;
}

sub parameter {
  my ($self, $key, $value) = @_;
  $self->set_parameter($key, $value) if $value;
  return $self->get_parameter($key);
}

sub scalex {
  my ($self, $scalex) = @_;
  
  if ($scalex) {
    $self->{'_scalex'} = $scalex;
    $self->texthelper->scalex($scalex);
  }
  
  return $self->{'_scalex'};
}

sub cache {
  my $self = shift;
  my $key  = shift;
  $self->{'_cache'}{$key} = shift if @_;
  return $self->{'_cache'}{$key}
}

# Delete all tracks where menu = no, and parent nodes if they are now empty
sub remove_disabled_menus {
  my ($self, $node) = @_;
  
  if (!$node) {
    $_->remove for grep $_->get('menu') eq 'no', $self->tree->leaves;
    $self->remove_disabled_menus($_) for $self->tree->nodes;
    return;
  }
  
  if ($node->get('node_type') !~ /^(track|option)$/ && !$node->has_child_nodes) {
    my $parent = $node->parent_node;
    $node->remove;
    $self->remove_disabled_menus($parent) if $parent && !scalar @{$parent->child_nodes};
  }
}

# create_menus - takes an array to configure the menus to be seen on the display
sub create_menus {
  my $self = shift;
  my $tree = $self->tree;
  
  foreach (@_) {
    my $menu = $self->menus->{$_};
    
    if (ref $menu) {
      my $parent = $tree->get_node($menu->[1]) || $tree->append_child($self->create_submenu($menu->[1], $self->menus->{$menu->[1]}));
      $parent->append_child($self->create_submenu($_, $menu->[0]));
    } else {
      $tree->append_child($self->create_submenu($_, $menu));
    }
  }
}

sub create_submenu {
  my ($self, $code, $caption, $options) = @_;
  
  my $details = {
    caption    => $caption, 
    node_type  => 'menu',
    %{$options || {}}
  };
  
  return $self->tree->create_node($code, $details);
}

sub create_track {
  my ($self, $code, $caption, $options) = @_;
  
  my $details = { 
    name      => $caption,
    node_type => 'track',
    %{$options || {}}
  };
  
  $details->{'strand'}    ||= 'b';      # Make sure we have a strand setting
  $details->{'display'}   ||= 'normal'; # Show unless we explicitly say no
  $details->{'renderers'} ||= [qw(off Off normal On)];
  $details->{'colours'}   ||= $self->species_defs->colour($options->{'colourset'}) if exists $options->{'colourset'};
  $details->{'glyphset'}  ||= $code;
  $details->{'caption'}   ||= $caption;
  
  return $self->tree->create_node($code, $details);
}

sub add_track { shift->add_tracks(shift, \@_); }

sub add_tracks {
  my $self     = shift;
  my $menu_key = shift;
  my $menu     = $self->get_node($menu_key);

  return unless $menu;

  foreach my $row (@_) {
    my ($key, $caption, $glyphset, $params) = @$row;
    my $node = $self->get_node($key);

    next if $node && $node->get('node_type') eq 'track';

    $params->{'glyphset'} = $glyphset;
    $menu->append($self->create_track($key, $caption, $params));
  }
}

sub create_option {
  my ($self, $code, $caption, $values, $renderers, $display) = @_;
  
  $values    ||= { off => 0, normal => 1 };
  $renderers ||= [ 'off', 'Off', 'normal', 'On' ];
  
  return $self->tree->create_node($code, {
    node_type => 'option',
    caption   => $caption,
    name      => $caption,
    values    => $values,
    renderers => $renderers,
    display   => $display || 'normal'
  });
}

sub add_option {
  my $self = shift;
  my $menu = $self->get_node('options');
  
  return unless $menu;
  
  $menu->append($self->create_option(@_));
}

sub add_options {
  my $self = shift;
  my $menu = $self->get_node(ref $_[0] ? 'options' : shift);
  
  return unless $menu;
  
  $menu->append($self->create_option(@$_)) for @_;
}

sub get_option {
  my ($self, $code, $key) = @_;
  my $node = $self->get_node($code);
  return $node ? $node->get($key || 'values')->{$node->get('display')} : 0;
}

# Order submenus alphabetically by caption
sub alphabetise_tracks {
  my ($self, $track, $menu, $key) = @_;
  $key ||= 'caption';
  
  my $name = $track->data->{$key};
  my ($after, $node_name);
  
  foreach (@{$menu->child_nodes}) {
    $node_name = $_->data->{$key};
    $after     = $_ if $node_name && $node_name lt $name;
  }
  
  if ($after) {
    $after->after($track);
  } else {
    $menu->prepend_child($track);
  }
}

sub load_user_tracks {
  my $self = shift;
  my $menu = $self->get_node('user_data');
  
  return unless $menu;
  
  my $hub      = $self->hub;
  my $session  = $hub->session;
  my $user     = $hub->user;
  my $das      = $hub->get_all_das;
  my $datahubs = $self->get_parameter('datahubs') == 1;
  my (%url_sources, %upload_sources);
  
  $self->_load_url_feature($menu);
  
  foreach my $source (sort { ($a->caption || $a->label) cmp ($b->caption || $b->label) } values %$das) {
    my $node = $self->get_node('das_' . $source->logic_name);

    next if     $node && $node->get('node_type') eq 'track';
    next unless $source->is_on($self->{'type'});
    
    $self->add_das_tracks('user_data', $source);
  }

  # Get the tracks that are temporarily stored - as "files" not in the DB....
  # Firstly "upload data" not yet committed to the database...
  # Then those attached as URLs to either the session or the User
  # Now we deal with the url sources... again flat file
  foreach my $entry ($session->get_data(type => 'url')) {
    next unless $entry->{'species'} eq $self->{'species'};
    
    $url_sources{"url_$entry->{'code'}"} = {
      source_type => 'session',
      source_name => $entry->{'name'} || $entry->{'url'},
      source_url  => $entry->{'url'},
      species     => $entry->{'species'},
      format      => $entry->{'format'},
      style       => $entry->{'style'},
      colour      => $entry->{'colour'},
      display     => $entry->{'display'},
      timestamp   => $entry->{'timestamp'} || time,
    };
  }
  
  foreach my $entry ($session->get_data(type => 'upload')) {
    next unless $entry->{'species'} eq $self->{'species'};
   
    if ($entry->{'analyses'}) {
      foreach my $analysis (split /, /, $entry->{'analyses'}) {
        $upload_sources{$analysis} = {
          source_name => $entry->{'name'},
          source_type => 'session',
          assembly    => $entry->{'assembly'},
          style       => $entry->{'style'},
        };
        
        $self->_compare_assemblies($entry, $session);
      }
    } elsif ($entry->{'species'} eq $self->{'species'} && !$entry->{'nonpositional'}) {
      my ($strand, $renderers) = $self->_user_track_settings($entry->{'style'}, $entry->{'format'});
      $strand = $entry->{'strand'} if $entry->{'strand'};
      
      $menu->append($self->create_track("upload_$entry->{'code'}", $entry->{'name'}, {
        external    => 'user',
        glyphset    => '_flat_file',
        colourset   => 'classes',
        sub_type    => 'tmp',
        file        => $entry->{'filename'},
        format      => $entry->{'format'},
        caption     => $entry->{'name'},
        renderers   => $renderers,
        description => 'Data that has been temporarily uploaded to the web server.',
        display     => 'off',
        strand      => $strand,
      }));
    }
  }
  
  if ($user) {
    my @groups = $user->get_groups;

    foreach my $entry (grep $_->species eq $self->{'species'}, $user->get_records('urls'), map $user->get_group_records($_, 'urls'), @groups) {
      $url_sources{'url_' . $entry->code} = {
        source_name => $entry->name || $entry->url,
        source_type => 'user', 
        source_url  => $entry->url,
        species     => $entry->species,
        format      => $entry->format,
        style       => $entry->style,
        colour      => $entry->colour,
        display     => 'off',
        timestamp   => $entry->timestamp,
      };
    }
    
    foreach my $entry (grep $_->species eq $self->{'species'}, $user->get_records('uploads'), map $user->get_group_records($_, 'uploads'), @groups) {
      my ($name, $assembly) = ($entry->name, $entry->assembly);
      
      foreach my $analysis (split /, /, $entry->analyses) {
        $upload_sources{$analysis} = {
          source_name => $name,
          source_type => 'user',
          assembly    => $assembly,
          style       => $entry->style,
        };
        
        $self->_compare_assemblies($entry, $session);
      }
    }
  }
  
  foreach my $code (sort { $url_sources{$a}{'source_name'} cmp $url_sources{$b}{'source_name'} } keys %url_sources) {
    my $add_method = lc "_add_$url_sources{$code}{'format'}_track";
    
    if ($self->can($add_method)) {
      $self->$add_method(
        key      => $code,
        menu     => $menu,
        source   => $url_sources{$code},
        external => 'user'
      );
    } elsif (lc $url_sources{$code}{'format'} eq 'datahub') {
      $self->_add_datahub($url_sources{$code}{'source_name'}, $url_sources{$code}{'source_url'}) if $datahubs;
    } else {
      $self->_add_flat_file_track($menu, 'url', $code, $url_sources{$code}{'source_name'},
        sprintf('
          Data retrieved from an external webserver. This data is attached to the %s, and comes from URL: %s',
          encode_entities($url_sources{$code}{'source_type'}), encode_entities($url_sources{$code}{'source_url'})
        ),
        url      => $url_sources{$code}{'source_url'},
        format   => $url_sources{$code}{'format'},
        style    => $url_sources{$code}{'style'},
        external => 'user',
      );
    }
  }
  
  # We now need to get a userdata adaptor to get the analysis info
  if (keys %upload_sources) {
    my $dbs        = EnsEMBL::Web::DBSQL::DBConnection->new($self->{'species'});
    my $dba        = $dbs->get_DBAdaptor('userdata');
    my $an_adaptor = $dba->get_adaptor('Analysis');
    my @tracks;
    
    foreach my $logic_name (keys %upload_sources) {
      my $analysis = $an_adaptor->fetch_by_logic_name($logic_name);
      
      next unless $analysis;
   
      $analysis->web_data->{'style'} ||= $upload_sources{$logic_name}{'style'};
     
      my ($strand, $renderers) = $self->_user_track_settings($analysis->web_data->{'style'}, $analysis->program_version);
      my $source_name = encode_entities($upload_sources{$logic_name}{'source_name'});
      my $description = encode_entities($analysis->description) || "User data from dataset $source_name";
      my $caption     = encode_entities($analysis->display_label);
         $caption     = "$source_name: $caption" unless $caption eq $upload_sources{$logic_name}{'source_name'};
         $strand      = $upload_sources{$logic_name}{'strand'} if $upload_sources{$logic_name}{'strand'};
      
      push @tracks, [ $logic_name, $caption, {
        external    => 'user',
        glyphset    => '_user_data',
        colourset   => 'classes',
        sub_type    => $upload_sources{$logic_name}{'source_type'} eq 'user' ? 'user' : 'tmp',
        renderers   => $renderers,
        source_name => $source_name,
        logic_name  => $logic_name,
        caption     => $caption,
        data_type   => $analysis->module,
        description => $description,
        display     => 'off',
        style       => $analysis->web_data,
        format      => $analysis->program_version,
        strand      => $strand,
      }];
    }
   
    $menu->append($self->create_track(@$_)) for sort { lc $a->[2]{'source_name'} cmp lc $b->[2]{'source_name'} || lc $a->[1] cmp lc $b->[1] } @tracks;
  }
 
  $ENV{'CACHE_TAGS'}{'user_data'} = sprintf 'USER_DATA[%s]', md5_hex(join '|', map $_->id, $menu->nodes) if $menu->has_child_nodes;
}

sub _add_datahub {
  my ($self, $menu_name, $url) = @_;
  my $parser   = Bio::EnsEMBL::ExternalData::DataHub::SourceParser->new({ timeout => 10, proxy => $self->hub->species_defs->ENSEMBL_WWW_PROXY });
  my $menu     = $self->tree->append_child($self->create_submenu($menu_name, $menu_name, { external => 1, datahub_menu => 1 }));
  my $hub_info = $parser->get_hub_info($url); ## Do we have data for this species?
  
  if ($hub_info->{'error'}) {
    warn "!!! COULD NOT CONTACT DATAHUB $url: $hub_info->{'error'}";
  } else {
    my $golden_path = $self->hub->species_defs->get_config($self->species, 'UCSC_GOLDEN_PATH');
    my $source_list = $hub_info->{'genomes'}{$golden_path} || [];
    
    return unless scalar @$source_list;
    
    ## Get tracks from hub
    my $datahub = $parser->parse($source_list);
    
    foreach my $node (@{$datahub->child_nodes}) {
      my $data = $node->data;
      
      if ($data->{'error'}) {
        warn "!!! COULD NOT PARSE CONFIG $data->{'file'}: $data->{'error'}";
      } elsif (!$node->has_child_nodes) {
        # No inheritance structure - assumes that the top level in the hub contains only tracks
        $self->_add_datahub_node($node->parent_node, $menu, $menu_name);
        last;
      } else {
        $self->_add_datahub_node($node, $menu, $menu_name);
      }
    }
  }
}

sub _add_datahub_node {
  my ($self, $node, $menu, $name) = @_;
  
  if ($node->has_child_nodes && $node->first_child->has_child_nodes) {
    $self->_add_datahub_node($_, $menu, $name) for @{$node->child_nodes};
  } else {
    my $config = {};
    my $n      = $node;
    
    do {
      my $data = $n->data;
      $config->{$_} ||= $data->{$_} for keys %$data;
    } while $n = $n->parent_node;
    
    $self->_add_datahub_tracks($node, $config, $menu, $name);
  }
}

sub _add_datahub_tracks {
  my ($self, $parent, $config, $menu, $name) = @_;
  my $hub    = $self->hub;
  my $data   = $parent->data;
  my $matrix = $config->{'dimensions'}{'x'} && $config->{'dimensions'}{'y'};
  my $link   = $config->{'description_url'} ? qq{<br /><a href="$config->{'description_url'}" rel="external">Go to track description on datahub</a>} : '';
  my $info   = $config->{'longLabel'} . $link;
  my %tracks;
  
  my %options = (
    menu_key     => $name,
    menu_name    => $name,
    submenu_key  => $self->tree->clean_id("${name}_$data->{'track'}", '\W'),
    submenu_name => $data->{'shortLabel'},
    datahub      => 1,
  );
  
  if ($matrix) {
    $options{'matrix_url'} = $hub->url('Config', { action => 'Matrix', function => $hub->action, partial => 1, menu => $options{'submenu_key'} });
    
    foreach my $subgroup (keys %$config) {
      next unless $subgroup =~ /subGroup\d/;
      
      foreach (qw(x y)) {
        if ($config->{$subgroup}{'name'} eq $config->{'dimensions'}{$_}) {
          $options{'axis_labels'}{$_} = { %{$config->{$subgroup}} }; # Make a deep copy so that the regex below doesn't affect the subgroup config
          s/_/ /g for values %{$options{'axis_labels'}{$_}};
        }
      }
      
      last if scalar keys %{$options{'axis_labels'}} == 2;
    }
    
    $options{'axes'} = { map { $_ => $options{'axis_labels'}{$_}{'label'} } qw(x y) };
  }
  
  my $submenu = $self->create_submenu($options{'submenu_key'}, $options{'submenu_name'}, {
    external => 1,
    ($matrix ? (
      menu   => 'matrix',
      url    => $options{'matrix_url'},
      matrix => {
        section     => $menu->data->{'caption'},
        header      => $options{'submenu_name'},
        description => $info,
        axes        => $options{'axes'},
      }
    ) : ())
  });
  
  $self->alphabetise_tracks($submenu, $menu);
  
  foreach (@{$parent->child_nodes}) {
    my $track        = $_->data;
    my $type         = ref $track->{'type'} eq 'HASH' ? uc $track->{'type'}{'format'} : uc $track->{'type'};
    my $squish       = $track->{'visibility'} eq 'squish' || $config->{'visibility'} eq 'squish'; # FIXME: make it inherit correctly
    (my $source_name = $track->{'shortLabel'}) =~ s/_/ /g;
    my $source       = {
      name        => $track->{'track'},
      source_name => $source_name,
      description => $track->{'longLabel'} . $link,
      source_url  => $track->{'bigDataUrl'},
      colour      => exists $track->{'color'} ? $track->{'color'} : undef,
      no_titles   => $type eq 'BIGWIG', # To improve browser speed don't display a zmenu for bigwigs
      squish      => $squish,
      %options
    };
    
    # Graph range - Track Hub default is 0-127
    if (exists $track->{'viewLimits'}) {
      $source->{'viewLimits'} = $track->{'viewLimits'};
    } elsif ($track->{'autoScale'} eq 'off') {
      $source->{'viewLimits'} = '0:127';
    }

    if (exists $track->{'maxHeightPixels'}) {
      $source->{'maxHeightPixels'} = $track->{'maxHeightPixels'};
    } elsif ($type eq 'BIGWIG' || $type eq 'BIGBED') {
      $source->{'maxHeightPixels'} = '64:32:16';
    }
    
    if ($matrix) {
      # Should really be shortLabel, but Encode is much better using longLabel (duplicate names using shortLabel)
      # The problem is that UCSC browser has a grouped set of tracks with a header above them. This
      # means the short label can be very non specific, because the header gives context of what type of
      # track it is. For Ensembl we need to have all the information in the track name / caption
      ($source->{'source_name'} = $track->{'longLabel'}) =~ s/_/ /g;
     
      $source->{'matrix'} = {
        menu   => $options{'submenu_key'},
        column => $options{'axis_labels'}{'x'}{$track->{'subGroups'}{$config->{'dimensions'}{'x'}}},
        row    => $options{'axis_labels'}{'y'}{$track->{'subGroups'}{$config->{'dimensions'}{'y'}}},
      };
      
      $source->{'column_data'} = { description => $info, no_subtrack_description => 1 };
    }
    
    $tracks{$type}{$source->{'name'}} = $source;
  }
  
  $self->load_file_format(lc, $tracks{$_}) for keys %tracks;
}

sub _add_datahub_extras_options {
  my ($self, %args) = @_;
  
  if (exists $args{'menu'}{'maxHeightPixels'} || exists $args{'source'}{'maxHeightPixels'}) {
    $args{'options'}{'maxHeightPixels'} = $args{'menu'}{'maxHeightPixels'} || $args{'source'}{'maxHeightPixels'};

    (my $default_height = $args{'options'}{'maxHeightPixels'}) =~ s/^.*:([0-9]*):.*$/$1/;
    
    $args{'options'}{'height'} = $default_height if $default_height > 0;
  }
  
  # Alternative rendering order for genome segmentation and similar
  if ($args{'source'}{'squish'}) {
    $args{'renderers'} = [
      'off',     'Off',
      'compact', 'Continuous',
      'normal',  'Separate',
      'labels',  'Separate with labels',
    ];
  }
  
  $args{'options'}{'viewLimits'} = $args{'menu'}{'viewLimits'} || $args{'source'}{'viewLimits'} if exists $args{'menu'}{'viewLimits'} || exists $args{'source'}{'viewLimits'};
  $args{'options'}{'no_titles'}  = $args{'menu'}{'no_titles'}  || $args{'source'}{'no_titles'}  if exists $args{'menu'}{'no_titles'}  || exists $args{'source'}{'no_titles'};
  $args{'options'}{'set'}        = $args{'source'}{'submenu_key'};
  $args{'options'}{'subset'}     = $self->tree->clean_id($args{'source'}{'submenu_key'}, '\W') unless $args{'source'}{'matrix'};
  $args{'options'}{$_}           = $args{'source'}{$_} for qw(datahub matrix column_data colour description);
  
  return %args;
}

sub _load_url_feature {
## Creates a temporary track based on a single line of a row-based data file,
## such as VEP output, e.g. 21 9678256 9678256 T/G 1 
  my ($self, $menu) = @_;
  return unless $menu;

  my $session_data = $self->hub->session->get_data('code' => 'custom_feature', type => 'custom');
  my ($format, $data);

  if ($self->hub->param('custom_feature')) {
    $format  = $self->hub->param('format');
    $data    = decode_entities($self->hub->param('custom_feature'));
    $session_data = {'code' => 'custom_feature', 'type' => 'custom', 
                      'format' => $format, 'data' => $data}; 
    $self->hub->session->set_data(%$session_data);
  }
  elsif ($session_data && ref($session_data) eq 'HASH' && $session_data->{'data'}) {
    $format = $session_data->{'format'};
    $data = $session_data->{'data'};
  }
  return unless ($data && $format);

  my ($strand, $renderers) = $self->_user_track_settings(undef, $format);
  my $file_info = $self->hub->species_defs->DATA_FORMAT_INFO;

  my $track = $self->create_track('custom_feature', 'Single feature', {
        external    => 'user',
        glyphset    => '_flat_file',
        colourset   => 'classes',
        sub_type    => 'single_feature',
        format      => $format,
        caption     => 'Single '.$file_info->{$format}{'label'}.' feature',
        renderers   => $renderers,
        description => 'A single feature that has been loaded via a hyperlink',
        display     => 'off',
        strand      => $strand,
        data        => $data,
  });
  $menu->append($track) if $track;
}

sub load_configured_bam    { shift->load_file_format('bam');    }
sub load_configured_bigbed { shift->load_file_format('bigbed'); }
sub load_configured_bigwig { shift->load_file_format('bigwig'); }
sub load_configured_vcf    { shift->load_file_format('vcf');    }

sub load_file_format {
  my ($self, $format, $sources) = @_;
  my $function = "_add_${format}_track";
  
  return unless $self->can($function);
  
  my $internal = !defined $sources;
     $sources  = $self->sd_call(sprintf 'ENSEMBL_INTERNAL_%s_SOURCES', uc $format) || {} unless defined $sources; # get the internal sources from config
  
  foreach my $source_name (sort keys %$sources) {
    # get the target menu 
    my $menu = $self->get_node($sources->{$source_name});
    my ($source, $view);
    
    if ($menu) {
      $source = $self->sd_call($source_name);
      $view   = $source->{'view'};
    } else {
      ## Probably an external datahub source
         $source       = $sources->{$source_name};
         $view         = $source->{'view'};
      my $menu_key     = $source->{'menu_key'};
      my $menu_name    = $source->{'menu_name'};
      my $submenu_key  = $source->{'submenu_key'};
      my $submenu_name = $source->{'submenu_name'};
      my $main_menu    = $self->get_node($menu_key) || $self->tree->append_child($self->create_submenu($menu_key, $menu_name, { external => 1, datahub_menu => !!$source->{'datahub'} }));
         $menu         = $self->get_node($submenu_key);
      
      if (!$menu) {
        $menu = $self->create_submenu($submenu_key, $submenu_name, { external => 1, ($source->{'matrix_url'} ? (menu => 'matrix', url => $source->{'matrix_url'}) : ()) });
        $self->alphabetise_tracks($menu, $main_menu);
      }
    }
    
    $self->$function(key => $source_name, menu => $menu, source => $source, description => $source->{'description'}, internal => $internal, view => $view) if $source;
  }
}

sub _add_bam_track {
  my $self = shift;
  my $desc = '
    The read end bars indicate the direction of the read and the colour indicates the type of read pair:
    Green = both mates are part of a proper pair; Blue = either this read is not paired, or its mate was not mapped; Red = this read is not properly paired.
  ';
  
  $self->_add_file_format_track(
    format      => 'BAM',
    description => $desc,
    renderers   => [
      'off',       'Off', 
      'normal',    'Normal', 
      'unlimited', 'Unlimited', 
      'histogram', 'Coverage only'
    ], 
    options => {
      external => 'external',
      sub_type => 'bam'
    },
    @_
  );
}

sub _add_bigbed_track {
  my ($self, %args) = @_;
 
  my $renderers = $args{'source'}{'renderers'} || [
    'off',     'Off', 
    'normal',  'Normal', 
    'labels',  'Labels',
    'compact', 'Compact',
  ];

  my $options = {
    external     => 'external',
    sub_type     => 'url',
    colourset    => 'feature',
    style        => $args{'source'}{'style'},
    addhiddenbgd => 1,
  };

  if ($args{'view'} && $args{'view'} =~ /peaks/i) {
    $options->{'join'} = 'off';  
  } else {
    push @$renderers, ('tiling', 'Wiggle plot');
  } 
  
  $self->_add_file_format_track(
    format      => 'BigBed',
    description => 'Bigbed file',
    renderers   => $renderers,
    options     => $options,
    %args,
  );
}

sub _add_bigwig_track {
  my ($self, %args) = @_;
  
  my $renderers = $args{'source'}{'renderers'} || [
    'off',     'Off',
    'tiling',  'Wiggle plot',
    'compact', 'Compact',
  ];
 
  my $options = {
    external     => 'external',
    sub_type     => 'bigwig',
    colour       => $args{'menu'}{'colour'} || $args{'source'}{'colour'} || 'red',
    addhiddenbgd => 1,
  };
  
  $self->_add_file_format_track(
    format    => 'BigWig', 
    renderers =>  $renderers,
    options   => $options,
    %args
  );
}

sub _add_vcf_track {
  shift->_add_file_format_track(
    format    => 'VCF', 
    renderers => [
      'off',       'Off',
      'histogram', 'Histogram',
      'compact',   'Compact'
    ], 
    options => {
      external   => 'external',
      sources    => undef,
      depth      => 0.5,
      bump_width => 0,
      colourset  => 'variation'
    },
    @_
  );
}

sub _add_flat_file_track {
  my ($self, $menu, $sub_type, $key, $name, $description, %options) = @_;
  
  $menu ||= $self->get_node('user_data');
  
  return unless $menu;
 
  my ($strand, $renderers) = $self->_user_track_settings($options{'style'}, $options{'format'});
  
  my $track = $self->create_track($key, $name, {
    display     => 'off',
    strand      => $strand,
    external    => 'external',
    glyphset    => '_flat_file',
    colourset   => 'classes',
    caption     => $name,
    sub_type    => $sub_type,
    renderers   => $renderers,
    description => $description,
    %options
  });
  
  $menu->append($track) if $track;
}

sub _add_file_format_track {
  my ($self, %args) = @_;
  my $menu = $args{'menu'} || $self->get_node('user_data');
  
  return unless $menu;
  
  %args = $self->_add_datahub_extras_options(%args) if $args{'source'}{'datahub'};
  
  my $type    = lc $args{'format'};
  my $article = $args{'format'} =~ /^[aeiou]/ ? 'an' : 'a';
  my ($desc, $url);
  
  if ($args{'internal'}) {
    $desc = "Data served from a $args{'format'} file: $args{'description'}";
    $url = join '/', $self->hub->species_defs->DATAFILE_BASE_PATH, lc $self->hub->species, $self->hub->species_defs->ASSEMBLY_NAME, $args{'source'}{'dir'}, $args{'source'}{'file'};
    $args{'options'}{'external'} = undef;
  } else {
    $desc = sprintf(
      'Data retrieved from %s %s file on an external webserver. %s This data is attached to the %s, and comes from URL: %s',
      $article,
      $args{'format'},
      $args{'description'},
      encode_entities($args{'source'}{'source_type'}), 
      encode_entities($args{'source'}{'source_url'})
    );
  }
  
  $self->generic_add($menu, undef, $args{'key'}, {}, {
    display     => 'off',
    strand      => 'f',
    format      => $args{'format'},
    glyphset    => $type,
    colourset   => $type,
    renderers   => $args{'renderers'},
    name        => $args{'source'}{'source_name'},
    caption     => exists($args{'source'}{'caption'}) ? $args{'source'}{'caption'} : $args{'source'}{'source_name'},
    url         => $url || $args{'source'}{'source_url'},
    description => $desc,
    %{$args{'options'}}
  });
}

sub _user_track_settings {
  my ($self, $style, $format) = @_;
  my ($strand, @user_renderers);
      
  if ($style =~ /^(wiggle|WIG)$/) {
    $strand         = 'r';
    @user_renderers = ('off', 'Off', 'tiling', 'Wiggle plot');
  } else {
    $strand         = uc($format) eq 'VEP_INPUT' ? 'f' : 'b'; 
    @user_renderers = (@{$self->{'alignment_renderers'}}, 'difference', 'Differences');
  }
  
  return ($strand, \@user_renderers);
}

sub _compare_assemblies {
  my ($self, $entry, $session) = @_;

  if ($entry->{'assembly'} && $entry->{'assembly'} ne $self->sd_call('ASSEMBLY_NAME')) {
    $session->add_data(
      type     => 'message',
      code     => 'userdata_assembly_mismatch',
      message  => "Sorry, track $entry->{'name'} is on an old assembly ($entry->{'assembly'}) and cannot be shown",
      function => '_error'
    );
  }
}

sub update_from_input {
  my $self  = shift;
  my $input = $self->hub->input;
  
  return $self->reset if $input->param('reset');
  
  my $diff   = $input->param('image_config');
  my $reload = 0;
  
  if ($diff) {
    my $track_reorder = 0;
    
    $diff = from_json($diff);
    $self->update_track_renderer($_, $diff->{$_}->{'renderer'}) for grep exists $diff->{$_}->{'renderer'}, keys %$diff;
    
    $reload        = $self->altered;
    $track_reorder = $self->update_track_order($diff) if $diff->{'track_order'};
    $reload      ||= $track_reorder;
    $self->update_favourite_tracks($diff);
  } else {
    my %favourites;
    
    foreach my $p ($input->param) {
      my $val = $input->param($p);
      
      if ($val =~ /favourite_(on|off)/) {
        $favourites{$p} = { favourite => $1 eq 'on' ? 1 : 0 };
      } else {
        $self->update_track_renderer($p, $val);
      }
    }
    
    $reload = $self->altered;
    
    $self->update_favourite_tracks(\%favourites) if scalar keys %favourites;
  }
  
  return $reload;
}

sub update_from_url {
  ## Tracks added "manually" in the URL (e.g. via a link)
  
  my ($self, @values) = @_;
  my $hub     = $self->hub;
  my $session = $hub->session;
  my $species = $hub->species;
  
  foreach my $v (@values) {
    my $format = $hub->param('format');
    my ($key, $renderer);
    
    if (uc $format eq 'DATAHUB') {
      $key = $v;
    } else {
      my @split = split /=/, $v;
      
      if (scalar @split > 1) {
        $renderer = pop @split;
        $key      = join '=', @split;
      } else {
        $key      = $split[0];
        $renderer = 'normal';
      }
    }
    
    if ($key =~ /^(\w+)[\.:](.*)$/) {
      my ($type, $p) = ($1, $2);
      
      if ($type eq 'url') {
        my $menu_name   = $hub->param('menu');
        my $all_formats = $hub->species_defs->DATA_FORMAT_INFO;
        
        if (!$format) {
          $p = uri_unescape($p);
          
          my @path = split(/\./, $p);
          my $ext  = $path[-1] eq 'gz' ? $path[-2] : $path[-1];
          
          while (my ($name, $info) = each %$all_formats) {
            if ($ext =~ /^$name$/i) {
              $format = $name;
              last;
            }  
          }
          if (!$format) {
            # Didn't match format name - now try checking format extensions
            while (my ($name, $info) = each %$all_formats) {
              if ($ext eq $info->{'ext'}) {
                $format = $name;
                last;
              }  
            }
          }
        }

        my $style = $all_formats->{lc $format}{'display'} eq 'graph' ? 'wiggle' : $format;
        my $code  = join '_', md5_hex("$species:$p"), $session->session_id;
        my $n;
        
        if ($menu_name) {
          $n = $menu_name;
        } else {
          $n = $p =~ /\/([^\/]+)\/*$/ ? $1 : 'un-named';
        }
        
        # Don't add if the URL or menu are the same as an existing track
        if ($session->get_data(type => 'url', code => $code)) {
          $session->add_data(
            type     => 'message',
            function => '_warning',
            code     => "duplicate_url_track_$code",
            message  => "You have already attached the URL $p. No changes have been made for this data source.",
          );
          
          next;
        } elsif (grep $_->{'name'} eq $n, $session->get_data(type => 'url')) {
          $session->add_data(
            type     => 'message',
            function => '_error',
            code     => "duplicate_url_track_$n",
            message  => qq{Sorry, the menu "$n" is already in use. Please change the value of "menu" in your URL and try again.},
          );
          
          next;
        }
        
        # We have to create a URL upload entry in the session
        $session->set_data(
          type    => 'url',
          url     => $p,
          species => $species,
          code    => $code, 
          name    => $n,
          format  => $format,
          style   => $style,
        );
        
        $session->add_data(
          type     => 'message',
          function => '_info',
          code     => 'url_data:' . md5_hex($p),
          message  => sprintf('Data has been attached to your display from the following URL: %s', encode_entities($p))
        );
        
        # We then have to create a node in the user_config
        if (uc $format eq 'DATAHUB') {
          $self->_add_datahub($n, $p);
        } else {
          $self->_add_flat_file_track(undef, 'url', "url_$code", $n, 
            sprintf('Data retrieved from an external webserver. This data is attached to the %s, and comes from URL: %s', encode_entities($n), encode_entities($p)),
            url   => $p,
            style => $style
          );

          $self->update_track_renderer("url_$code", $renderer);
        }       
      } elsif ($type eq 'das') {
        $p = uri_unescape($p);

        my $logic_name = $session->add_das_from_string($p, $self->{'type'}, { display => $renderer });

        if ($logic_name) {
          $session->add_data(
            type     => 'message',
            function => '_info',
            code     => 'das:' . md5_hex($p),
            message  => sprintf('You have attached a DAS source with DSN: %s %s.', encode_entities($p), $self->get_node("das_$logic_name") ? 'to this display' : 'but it cannot be displayed on the specified image')
          );
        }
      }
    } else {
      $self->update_track_renderer($key, $renderer, $hub->param('toggle_tracks'));
    }
  }
  
  if ($self->altered) {
    $session->add_data(
      type     => 'message',
      function => '_info',
      code     => 'image_config',
      message  => 'The link you followed has made changes to the tracks displayed on this page.',
    );
  }
}

sub update_track_renderer {
  my ($self, $key, $renderer, $on_off) = @_;
  my $node = $self->get_node($key);
  
  return unless $node;
  
  my $renderers = $node->data->{'renderers'};
  
  return unless $renderers;
  
  my %valid = @$renderers;
  my $flag  = 0;

  ## Set renderer to something sensible if user has specified invalid one. 'off' is usually first option, so take next one
  $renderer = $valid{'normal'} ? 'normal' : $renderers->[2] if $renderer ne 'off' && !$valid{$renderer};

  # if $on_off == 1, only allow track enabling/disabling. Don't allow enabled tracks' renderer to be changed.
  $flag += $node->set_user('display', $renderer) if (!$on_off || $renderer eq 'off' || $node->get('display') eq 'off');

  $self->altered = 1 if $flag;
}

sub update_favourite_tracks {
  my ($self, $diff) = @_;
  my $hub     = $self->hub;
  my $session = $hub->session;
  my $user    = $hub->user;
  my %args    = ( type => 'favourite_tracks', code => 'favourite_tracks' );
  
  $args{'tracks'} = $self->get_favourite_tracks;
  
  foreach (grep exists $diff->{$_}->{'favourite'}, keys %$diff) {
    if ($diff->{$_}->{'favourite'} == 1) {
      $args{'tracks'}{$_} = 1;
    } elsif (exists $diff->{$_}->{'favourite'}) {
      delete $args{'tracks'}{$_};
    }
  }
  
  if (scalar keys %{$args{'tracks'}}) {
    $session->set_data(%args);
  } else {
    delete $args{'tracks'};
    $session->purge_data(%args);
  }
  
  $user->set_favourite_tracks($args{'tracks'}) if $user;
  
  delete $self->{'favourite_tracks'};
}

sub update_track_order {
  my ($self, $diff) = @_;
  my $species    = $self->species;
  my $node       = $self->get_node('track_order');
  $self->altered = $node->set_user($species, { %{$node->get($species) || {}}, %{$diff->{'track_order'}} });
  return $self->altered if $self->get_parameter('sortable_tracks') ne 'drag';
}

sub reset {
  my $self  = shift;
  my $reset = $self->hub->input->param('reset');
  my ($tracks, $order) = $reset eq 'all' ? (1, 1) : $reset eq 'track_order' ? (0, 1) : (1, 0);
  
  if ($tracks) {
    my $tree = $self->tree;
    
    foreach my $node ($tree, $tree->nodes) {
      my $user_data = $node->{'user_data'};
      
      foreach (keys %$user_data) {
        $self->altered = 1 if $user_data->{$_}{'display'};
        delete $user_data->{$_}{'display'};
        delete $user_data->{$_} unless scalar keys %{$user_data->{$_}};
      }
    }
  }
  
  if ($order) {
    my $node    = $self->get_node('track_order');
    my $species = $self->species;
    
    if ($node->{'user_data'}{'track_order'}{$species}) {
      $self->altered = 1;
      delete $node->{'user_data'}{'track_order'}{$species};
      delete $node->{'user_data'}{'track_order'} unless scalar keys %{$node->{'user_data'}{'track_order'}};
    }
  }
}

sub get_track_key {
  my ($self, $prefix, $obj) = @_;
  
  return if $obj->gene && $obj->gene->isa('Bio::EnsEMBL::ArchiveStableId');

  my $logic_name = $obj->gene ? $obj->gene->analysis->logic_name : $obj->analysis->logic_name;
  my $db         = $obj->get_db;
  my $db_key     = 'DATABASE_' . uc $db;
  my $key        = $self->databases->{$db_key}{'tables'}{'gene'}{'analyses'}{lc $logic_name}{'web'}{'key'} || lc $logic_name;
  return join '_', $prefix, $db, $key;
}

sub modify_configs {
  my ($self, $nodes, $config) = @_;
  
  foreach my $node (map { $self->get_node($_) || () } @$nodes) {
    foreach my $n (grep $_->get('node_type') eq 'track', $node, $node->nodes) {
      $n->set($_, $config->{$_}) for keys %$config;
    }
  }
}

sub _update_missing {
  my ($self, $object) = @_;
  my $species_defs    = $self->species_defs;
  my $count_missing   = grep { $_->get('display') eq 'off' || !$_->get('display') } $self->get_tracks; 
  my $missing         = $self->get_node('missing');
  $missing->set('text', $count_missing > 0 ? "There are currently $count_missing tracks turned off." : 'All tracks are turned on') if $missing;
  
  my $info = sprintf(
    '%s %s version %s.%s (%s) %s: %s - %s',
    $species_defs->ENSEMBL_SITETYPE,
    $species_defs->SPECIES_BIO_NAME,
    $species_defs->ENSEMBL_VERSION,
    $species_defs->SPECIES_RELEASE_VERSION,
    $species_defs->ASSEMBLY_NAME,
    $object->seq_region_type_and_name,
    $object->thousandify($object->seq_region_start),
    $object->thousandify($object->seq_region_end)
  );

  my $information = $self->get_node('info');
  $information->set('text', $info) if $information;
  
  return { count => $count_missing, information => $info };
}

# load_tracks - loads in various database derived tracks; 
# loop through core like dbs, compara like dbs, funcgen like dbs, variation like dbs
sub load_tracks { 
  my $self         = shift;
  my $species      = $self->{'species'};
  my $species_defs = $self->species_defs;
  my $dbs_hash     = $self->databases;

  my %data_types = (
    core => [
      'add_dna_align_features',     # Add to cDNA/mRNA, est, RNA, other_alignment trees
      'add_data_files',             # Add to gene/rnaseq tree
#     'add_ditag_features',         # Add to ditag_feature tree
      'add_genes',                  # Add to gene, transcript, align_slice_transcript, tsv_transcript trees
      'add_trans_associated',       # Add to features associated with transcripts
      'add_marker_features',        # Add to marker tree
      'add_qtl_features',           # Add to marker tree
      'add_genome_attribs',         # Add to genome_attribs tree
      'add_misc_features',          # Add to misc_feature tree
      'add_prediction_transcripts', # Add to prediction_transcript tree
      'add_protein_align_features', # Add to protein_align_feature_tree
      'add_protein_features',       # Add to protein_feature_tree
      'add_repeat_features',        # Add to repeat_feature tree
      'add_simple_features',        # Add to simple_feature tree
      'add_decorations'
    ],
    compara => [
      'add_synteny',                # Add to synteny tree
      'add_alignments'              # Add to compara_align tree
    ],
    funcgen => [
      'add_regulation_builds',      # Add to regulation_feature tree
      'add_regulation_features',    # Add to regulation_feature tree
      'add_oligo_probes'            # Add to oligo tree
    ],
    variation => [
      'add_sequence_variations',          # Add to variation_feature tree
      'add_structural_variations',        # Add to variation_feature tree
      'add_copy_number_variant_probes',   # Add to variation_feature tree
      'add_phenotypes',                   # Add to variation_feature tree
      'add_recombination',                # Moves recombination menu to the end of the variation_feature tree
      'add_somatic_mutations',            # Add to somatic tree
      'add_somatic_structural_variations' # Add to somatic tree
    ],
  );
  
  foreach my $type (keys %data_types) {
    my ($check, $databases) = $type eq 'compara' ? ($species_defs->multi_hash, $species_defs->compara_like_databases) : ($dbs_hash, $self->sd_call("${type}_like_databases"));
    
    foreach my $db (grep exists $check->{$_}, @{$databases || []}) {
      my $key = lc substr $db, 9;
      $self->$_($key, $check->{$db}{'tables'} || $check->{$db}, $species) for @{$data_types{$type}}; # Look through tables in databases and add data from each one      
    }
  }
  
  $self->add_options('information', [ 'opt_empty_tracks', 'Display empty tracks', undef, undef, 'off' ]) unless $self->get_parameter('opt_empty_tracks') eq '0';
  $self->add_options('information', [ 'opt_highlight_feature', 'Highlight current feature', undef, undef, 'normal' ]);
  $self->tree->append_child($self->create_option('track_order')) if $self->get_parameter('sortable_tracks');
}

sub load_configured_das {
  my $self          = shift;
  my $extra         = ref $_[0] eq 'HASH' ? shift : {};
  my %allowed_menus = map { $_ => 1 } @_;
  my $all_menus     = !scalar @_;
  my @adding;
  my %seen;
  
  foreach my $source (sort { $a->caption cmp $b->caption } values %{$self->species_defs->get_all_das}) {
    next unless $source->is_on($self->{'type'});
    
    my ($category, $sub_category) = split ' ', $source->category;
    
    if ($category == 1) {
      $self->add_das_tracks('external_data', $source, $extra); # Unconfigured, will go into External data section
      next;
    }
    
    next unless $all_menus || $allowed_menus{$category};
    
    my $menu = $self->get_node($category);
    my $key;
    
    if (!$menu && grep { $category eq $_ } @{$self->{'transcript_types'}}) {
      foreach (@{$self->{'transcript_types'}}) {
        $category = $_ and last if $menu = $self->get_node($_);
      }
    }
    
    if (!$menu) {
      push @{$adding[0]}, $category unless $seen{$category}++;
      push @{$adding[1]{$category}}, $sub_category if $sub_category && !$seen{$sub_category}++;
      next;
    }
    
    if ($sub_category) {
      $key  = join '_', $category, lc $sub_category;
      
      my $sub_menu = $menu->get_node($key);
      
      if (!$sub_menu) {
        push @{$adding[1]{$category}}, $sub_category unless $seen{$sub_category}++;
        next;
      }
      
      if ($sub_menu && grep !$_->get('external'), @{$sub_menu->child_nodes}) {
        $menu = $sub_menu;
        $key  = "${key}_external";
      }
    } else {
      $key = "${category}_external";
    }
    
    $menu->append($self->create_submenu($key, 'External data', { external => 1 })) if $menu && !$menu->get_node($key);
    
    $self->add_das_tracks($key, $source, $extra);
  }
  
  # Add new menus, then run the function again - ensures that everything is printed in the right place
  if (scalar @adding) {
    my $external  = $self->get_node('external_data');
    my $menus     = $self->menus;
    my @new_menus = @{$adding[0] || []};
       %seen      = map { $_ => 1 } @new_menus;
    
    foreach (@new_menus) {
      my $parent = ref $menus->{$_} ? $self->get_node($menus->{$_}[1]) : undef;
      my $menu   = $self->get_node($_);
      
      $self->create_menus($_) unless $menu;
      
      $menu = $self->get_node($_);
      
      next unless $menu && $external;
      
      $external->after(ref $menus->{$_} ? $self->get_node($menus->{$_}[1]) : $menu) unless $parent;
    }
    
    foreach my $k (keys %{$adding[1]}) {
      $self->create_menus(@{$adding[1]{$k}});
      
      foreach (@{$adding[1]{$k}}) {
        my $key      = join '_', $k, lc $_;
        my $menu     = $self->get_node($k);
        my $sub_menu = $menu->get_node($key);
        
        if (!$sub_menu) {
          (my $caption = $_) =~ s/_/ /g;
          $menu->append($self->create_submenu($key, $caption));
        }
      }
      
      push @new_menus, $k unless $seen{$_};
    }
    
    $self->load_configured_das($extra, @new_menus) if scalar @new_menus;
  }
}

# Attach all das sources from an image config
sub attach_das {
  my $self      = shift;

  my @das_nodes = map { $_->get('glyphset') eq '_das' && $_->get('display') ne 'off' ? @{$_->get('logic_names')||[]} : () } $self->tree->nodes; # Look for all das sources which are configured and turned on

  return unless @das_nodes; # Return if no sources to be drawn
  
  my $hub         = $self->hub;
  my %T           = %{$hub->get_all_das}; # Check to see if they really exists, and get entries from get_all_das call
  my @das_sources = @T{@das_nodes};

  return unless @das_sources; # Return if no sources exist
  
  my $species_defs = $hub->species_defs;

  # Cache the DAS Coordinator object (with key das_coord)
  $self->cache('das_coord',  
    Bio::EnsEMBL::ExternalData::DAS::Coordinator->new(
      -sources => \@das_sources,
      -proxy   => $species_defs->ENSEMBL_WWW_PROXY,
      -noproxy => $species_defs->ENSEMBL_NO_PROXY,
      -timeout => $species_defs->ENSEMBL_DAS_TIMEOUT
    )
  );
}

sub _merge {
  my ($self, $_sub_tree, $sub_type) = @_;
  my $tree        = $_sub_tree->{'analyses'};
  my $config_name = $self->{'type'};
  my $data        = {};
  
  foreach my $analysis (keys %$tree){ 
    my $sub_tree = $tree->{$analysis}; 
    
    next unless $sub_tree->{'disp'}; # Don't include non-displayable tracks
    next if exists $sub_tree->{'web'}{$sub_type}{'do_not_display'};
    
    my $key = $sub_tree->{'web'}{'key'} || $analysis;
    
    foreach (grep $_ ne 'desc', keys %{$sub_tree->{'web'} || {}}) {
      if ($_ eq 'default') {
        $data->{$key}{'display'} ||= ref $sub_tree->{'web'}{$_} eq 'HASH' ? $sub_tree->{'web'}{$_}{$config_name} : $sub_tree->{'web'}{$_};
      } else {
        $data->{$key}{$_} ||= $sub_tree->{'web'}{$_}; # Longer form for help and configuration
      }
    }
    
    if ($sub_tree->{'web'}{'key'}) {
      if ($sub_tree->{'desc'}) {
        $data->{$key}{'multiple'}      = "This track comprises multiple analyses;" if $data->{$key}{'description'};
        $data->{$key}{'description'} ||= '';
        $data->{$key}{'description'}  .= ($data->{$key}{'description'} ? '; ' : '') . $sub_tree->{'desc'};
      }
    } else {
      $data->{$key}{'description'} = $sub_tree->{'desc'};
    }
    
    $data->{$key}{'format'} = $sub_tree->{'format'};
    
    push @{$data->{$key}{'logic_names'}}, $analysis;
  }
  
  foreach my $key (keys %$data) {
    $data->{$key}{'name'}      ||= $tree->{$key}{'name'};
    $data->{$key}{'caption'}   ||= $data->{$key}{'name'} || $tree->{$key}{'name'};
    $data->{$key}{'display'}   ||= 'off';
    $data->{$key}{'strand'}    ||= 'r';
    $data->{$key}{'description'} = "$data->{$key}{'multiple'} $data->{$key}{'description'}" if $data->{$key}{'multiple'};
  }
  
  return ([ sort { $data->{$a}{'name'} cmp $data->{$b}{'name'} } keys %$data ], $data);
}

sub generic_add {
  my ($self, $menu, $key, $name, $data, $options) = @_;
  
  $data = {
    %$data,
    db        => $key,
    renderers => [ 'off', 'Off', 'normal', 'On' ],
    %$options
  };
  
  $self->add_matrix($data, $menu) if $data->{'matrix'};
  
  return $menu->append($self->create_track($name, $data->{'name'}, $data));
}

sub add_matrix {
  my ($self, $data, $menu) = @_;
  my $menu_data    = $menu->data;
  my $matrix       = $data->{'matrix'};
  my $caption      = $data->{'caption'};
  my $after        = $data->{'track_after'};
  my $column       = $matrix->{'column'};
  my $subset       = $matrix->{'menu'};
  my @rows         = $matrix->{'rows'} ? @{$matrix->{'rows'}} : $matrix;
  my $column_key   = $self->tree->clean_id("${subset}_$column");
  my $column_track = $self->get_node($column_key);
  
  if (!($column_track && $column_track->parent_node)) {
    $column_track = $self->create_track($column_key, $data->{'track_name'} || $column, {
      renderers   => $data->{'renderers'},
      label_x     => $column,
      display     => 'off',
      subset      => $subset,
      track_after => $after,
      $matrix->{'row'} ? (matrix => 'column') : (),
      %{$data->{'column_data'} || {}}
    });
    
    $self->alphabetise_tracks($column_track, $menu, 'label_x');
  }
  
  if ($matrix->{'row'}) {
    push @{$column_track->data->{'subtrack_list'}}, [ $caption, $column_track->data->{'no_subtrack_description'} ? () : $data->{'description'} ];
    $data->{'option_key'} = $self->tree->clean_id("${subset}_${column}_$matrix->{'row'}");
  }
  
  $data->{'column_key'}  = $column_key;
  $data->{'menu'}        = 'matrix_subtrack';
  $data->{'source_name'} = $data->{'name'};
  $data->{'display'}     = 'default';
  
  if (!$menu_data->{'matrix'}) {
    my $hub = $self->hub;
    
    $menu_data->{'menu'}   = 'matrix';
    $menu_data->{'url'}    = $hub->url('Config', { action => 'Matrix', function => $hub->action, partial => 1, menu => $menu->id });
    $menu_data->{'matrix'} = {
      section => $menu->parent_node->data->{'caption'},
      header  => $menu_data->{'caption'},
    }
  }
  
  foreach (@rows) {
    my $option_key = $self->tree->clean_id("${subset}_${column}_$_->{'row'}");
    
    $column_track->append($self->create_track($option_key, $_->{'row'}, {
      node_type => 'option',
      menu      => 'no',
      display   => $_->{'on'} ? 'on' : 'off',
      renderers => [qw(on on off off)],
      caption   => "$column - $_->{'row'}",
    }));
    
    $menu_data->{'matrix'}{'rows'}{$_->{'row'}} ||= { id => $_->{'row'}, group => $_->{'group'}, group_order => $_->{'group_order'} };
  }
  
  return $column_track;
}

sub add_das_tracks {
  my ($self, $menu, $source, $extra) = @_;
  my $node = $self->get_node($menu); 
  
  if (!$node && grep { $menu eq "${_}_external" } @{$self->{'transcript_types'}}) {
    for (@{$self->{'transcript_types'}}) {
      $node = $self->get_node("${_}_external");
      last if $node;
    }
  }
  
  $node ||= $self->get_node('external_data'); 
  
  return unless $node;
  
  my $caption  = $source->caption || $source->label;
  my $desc     = $source->description;
  my $homepage = $source->homepage;
  
  $desc .= sprintf ' [<a href="%s" rel="external">Homepage</a>]', $homepage if $homepage;
  
  my $track = $self->create_track('das_' . $source->logic_name, $source->label, {
    %{$extra || {}},
    external    => 'external',
    glyphset    => '_das',
    display     => 'off',
    logic_names => [ $source->logic_name ],
    caption     => $caption,
    description => $desc,
    renderers   => [
      'off',      'Off', 
      'nolabels', 'No labels', 
      'normal',   'Normal', 
      'labels',   'Labels'
    ],
  });
  
  if ($track) {
    $node->append($track);
    $self->has_das ||= 1;
  }
}

#----------------------------------------------------------------------#
# Functions to add tracks from core like databases                     #
#----------------------------------------------------------------------#

# add_dna_align_features
# loop through all core databases - and attach the dna align
# features from the dna_align_feature tables...
# these are added to one of five menus: transcript, cdna/mrna, est, rna, other
# depending whats in the web_data column in the database
sub add_dna_align_features {
  my ($self, $key, $hashref) = @_;
  
  return unless $self->get_node('dna_align_cdna');
  
  my ($keys, $data) = $self->_merge($hashref->{'dna_align_feature'}, 'dna_align_feature');
  
  foreach my $key_2 (@$keys) {
    my $k    = $data->{$key_2}{'type'} || 'other';
    my $menu = ($k =~ /rnaseq|simple|transcript/) ? $self->tree->get_node($k) : $self->tree->get_node("dna_align_$k");
    
    if ($menu) {
      my $alignment_renderers = ['off','Off'];
      
      $alignment_renderers = [ @{$self->{'alignment_renderers'}} ] unless($data->{$key_2}{'no_default_renderers'});
            
      if (my @other_renderers = @{$data->{$key_2}{'additional_renderers'} || [] }) {
        my $i = 0;
        while ($i < scalar(@other_renderers)) {
          splice @$alignment_renderers, $i+2, 0, $other_renderers[$i];
          splice @$alignment_renderers, $i+3, 0, $other_renderers[$i+1];
          $i += 2;
        }
      }
      
      # my $display = (grep { $data->{$key_2}{'display'} eq $_ } @$alignment_renderers )             ? $data->{$key_2}{'display'}
      #             : (grep { $data->{$key_2}{'display'} eq $_ } @{$self->{'alignment_renderers'}} ) ? $data->{$key_2}{'display'}
      #             : 'off'; # needed because the same logic_name can be a gene and an alignment

      my $display  = $data->{$key_2}{'display'} ? $data->{$key_2}{'display'} : 'off';
      my $glyphset = '_alignment';
      my $strand   = 'b';
      
      if ($key_2 eq 'alt_seq_mapping') {
        $display             = 'simple';
        $alignment_renderers = [ 'off', 'Off', 'normal', 'On' ];  
        $glyphset            = 'patch_ref_alignment';
        $strand              = 'f';
      }
      
      $self->generic_add($menu, $key, "dna_align_${key}_$key_2", $data->{$key_2}, {
        glyphset  => $glyphset,
        sub_type  => lc $k,
        colourset => 'feature',
        display   => $display,
        renderers => $alignment_renderers,
        strand    => $strand,
      });
    }
  }
  
  $self->add_track('information', 'diff_legend', 'Alignment Difference Legend', 'diff_legend', { strand => 'r' });
}

sub add_data_files {
  my ($self, $key, $hashref) = @_;
  my $menu = $self->tree->get_node('rnaseq');

  return unless $menu;
  
  my ($keys, $data) = $self->_merge($hashref->{'data_file'});
  
  foreach (@$keys) {
    $self->generic_add($menu, $key, "data_file_${key}_$_", $data->{$_}, { 
      glyphset  => $data->{$_}{'format'}     || '_alignment', 
      colourset => $data->{$_}{'colour_key'} || 'feature',
      strand    => 'f',
      renderers => [
        'off',       'Off', 
        'normal',    'Normal', 
        'unlimited', 'Unlimited', 
        'histogram', 'Coverage only'
      ], 
    });
  }
}

sub add_ditag_features {
  my ($self, $key, $hashref) = @_;
  my $menu = $self->tree->get_node('ditag');
  
  return unless $menu;
  
  my ($keys, $data) = $self->_merge($hashref->{'ditag_feature'});
  $self->generic_add($menu, $key, "ditag_${key}_$_", $data->{$_}, { glyphset => '_ditag', strand => 'b' }) for @$keys;
}

# add_genes
# loop through all core databases - and attach the gene
# features from the gene tables...
# there are a number of menus sub-types these are added to:
# * transcript              # ordinary transcripts
# * alignslice_transcript   # transcripts in align slice co-ordinates
# * tse_transcript          # transcripts in collapsed intro co-ords
# * tsv_transcript          # transcripts in collapsed intro co-ords
# * gsv_transcript          # transcripts in collapsed gene co-ords
# depending on which menus are configured
sub add_genes {
  my ($self, $key, $hashref, $species) = @_;

  # Gene features end up in each of these menus
  return unless grep $self->get_node($_), @{$self->{'transcript_types'}};

  my ($keys, $data) = $self->_merge($hashref->{'gene'}, 'gene');
  my $colours       = $self->species_defs->colour('gene');
  my $flag          = 0;

  my $renderers = [
          'off',                     'Off',
          'gene_nolabel',            'No exon structure without labels',
          'gene_label',              'No exon structure with labels',
          'transcript_nolabel',      'Expanded without labels',
          'transcript_label',        'Expanded with labels',
          'collapsed_nolabel',       'Collapsed without labels',
          'collapsed_label',         'Collapsed with labels',
          'transcript_label_coding', 'Coding transcripts only (in coding genes)',          
        ];
        
  #push($renderers, 'transcript_gencode_basic','GENCODE basic'); # if($species eq "Homo_sapiens" || $species eq "Mus_musculus");  #only human and mouse have this renderer enable for now
     
  foreach my $type (@{$self->{'transcript_types'}}) {  
    my $menu = $self->get_node($type);
    next unless $menu;

    foreach my $key2 (@$keys) {
      my $t = $type;

      # force genes into a seperate menu if so specified in web_data (ie rna-seq); unless you're on a transcript page that is
      if ($data->{$key2}{'type'}){
        unless (ref($self) =~ /transcript/) {
          $t = $data->{$key2}{'type'};
        }
      }

      my $menu = $self->get_node($t);
      
      next unless $menu;
      $self->generic_add($menu, $key, "${t}_${key}_$key2", $data->{$key2}, {
        glyphset  => ($t =~ /_/ ? '' : '_') . $type, # QUICK HACK
        colours   => $colours,
        strand    => $t eq 'gene' ? 'r' : 'b',
        renderers => $t eq 'transcript' ? $renderers : $t eq 'rnaseq' ? [
         'off',                'Off',
         'transcript_nolabel', 'Expanded without labels',
         'transcript_label',   'Expanded with labels',
        ] : [
         'off',          'Off',
         'gene_nolabel', 'No labels', 
         'gene_label',   'With labels'
        ]
      });
      $flag = 1;
    }
  }

  # Need to add the gene menu track here
  $self->add_track('information', 'gene_legend', 'Gene Legend', 'gene_legend', { strand => 'r' }) if $flag;
}

sub add_trans_associated {
  my ($self, $key, $hashref) = @_;
  my $menu = $self->get_node('trans_associated');
  
  return unless $menu;
  
  my ($keys, $data) = $self->_merge($hashref->{'simple_feature'});
  $self->generic_add($menu, $key, "simple_${key}_$_", $data->{$_}, { glyphset => '_simple', colourset => 'simple' }) for grep $data->{$_}{'transcript_associated'}, @$keys;  
}

sub add_marker_features {
  my($self, $key, $hashref) = @_;
  my $menu = $self->get_node('marker');
  
  return unless $menu;
  
  my ($keys, $data) = $self->_merge($hashref->{'marker_feature'});
  my $colours = $self->species_defs->colour('marker');
  
  foreach (@$keys) {
    $self->generic_add($menu, $key, "marker_${key}_$_", $data->{$_}, {
      glyphset => '_marker',
      labels   => 'on',
      colours  => $colours,
      strand   => 'r',
    });
  }
}

sub add_qtl_features {
  my ($self, $key, $hashref) = @_;
  
  my $menu = $self->get_node('marker');
  
  return unless $menu && $hashref->{'qtl'} && $hashref->{'qtl'}{'rows'} > 0;
  
  $menu->append($self->create_track("qtl_$key", 'QTLs', {
    db          => $key,
    glyphset    => '_qtl',
    caption     => 'QTLs',
    colourset   => 'qtl',
    description => 'Quantative trait loci',
    display     => 'normal',
    renderers   => [ 'off', 'Off', 'normal', 'On' ],
    strand      => 'r',
  }));
}

sub add_genome_attribs {
  my ($self, $key, $hashref) = @_;
  my $menu = $self->get_node('genome_attribs');
  
  return unless $menu;
 
  my $default_tracks = {}; 
  my $config_name = $self->{'type'};
  my $data        = $hashref->{'genome_attribs'}{'sets'}; # Different loop - no analyses - just misc_sets
  
  foreach (sort { $data->{$a}{'name'} cmp $data->{$b}{'name'} } keys %$data) {
    next if $_ eq 'NoAnnotation' || $default_tracks->{$config_name}{$_}{'available'} eq 'no';
    
    $self->generic_add($menu, $key, "genome_attribs_${key}_$_", $data->{$_}, {
      glyphset          => '_clone',
      set               => $_,
      colourset         => 'clone',
      caption           => $data->{$_}{'name'},
      description       => $data->{$_}{'desc'},
      strand            => 'r',
      display           => $default_tracks->{$config_name}{$_}{'default'} || $data->{$_}{'display'} || 'off',
      outline_threshold => $default_tracks->{$config_name}{$_}{'threshold'} eq 'no' ? undef : 350000,
    });
  }
}

sub add_misc_features {
  my ($self, $key, $hashref) = @_;
  my $menu = $self->get_node('misc_feature');
  
  return unless $menu;
  
  # set some defaults and available tracks
  my $default_tracks = {
    cytoview   => {
      tilepath => { default   => 'normal' },
      encode   => { threshold => 'no'     }
    },
    contigviewbottom => {
      ntctgs => { available => 'no' },
      encode => { threshold => 'no' }
    }
  };
  
  my $config_name = $self->{'type'};
  my $data        = $hashref->{'misc_feature'}{'sets'}; # Different loop - no analyses - just misc_sets
  
  foreach (sort { $data->{$a}{'name'} cmp $data->{$b}{'name'} } keys %$data) {
    next if $_ eq 'NoAnnotation' || $default_tracks->{$config_name}{$_}{'available'} eq 'no';
    
    $self->generic_add($menu, $key, "misc_feature_${key}_$_", $data->{$_}, {
      glyphset          => '_clone',
      set               => $_,
      colourset         => 'clone',
      caption           => $data->{$_}{'name'},
      description       => $data->{$_}{'desc'},
      strand            => 'r',
      display           => $default_tracks->{$config_name}{$_}{'default'} || $data->{$_}{'display'} || 'off',
      outline_threshold => $default_tracks->{$config_name}{$_}{'threshold'} eq 'no' ? undef : 350000,
    });
  }
}

sub add_prediction_transcripts {
  my ($self, $key, $hashref) = @_;
  my $menu = $self->get_node('prediction');
  
  return unless $menu;
  
  my ($keys, $data) = $self->_merge($hashref->{'prediction_transcript'});
  
  foreach (@$keys) {
    $self->generic_add($menu, $key, "transcript_${key}_$_", $data->{$_}, {
      glyphset   => '_prediction_transcript',
      colourset  => 'prediction',
      label_key  => '[display_label]',
      colour_key => lc $_,
      renderers  => [ 'off', 'Off', 'transcript_nolabel', 'No labels', 'transcript_label', 'With labels' ],
      strand     => 'b',
    });
  }
}

# add_protein_align_features
# loop through all core databases - and attach the protein align
# features from the protein_align_feature tables...
sub add_protein_align_features {
  my ($self, $key, $hashref) = @_;
  my $menu = $self->tree->get_node('protein_align');
  
  return unless $menu;
  
  my ($keys, $data) = $self->_merge($hashref->{'protein_align_feature'}, 'protein_align_feature');
  
  foreach my $key_2 (@$keys) {
    # needed because the same logic_name can be a gene and an alignment, need to fix default rederer  the web_data
    my $display = (grep { $data->{$key_2}{'display'} eq $_ } @{$self->{'alignment_renderers'}}) ? $data->{$key_2}{'display'} : 'off';
    
    $self->generic_add($menu, $key, "protein_${key}_$key_2", $data->{$key_2}, {
      glyphset    => '_alignment',
      sub_type    => 'protein',
      colourset   => 'feature',
      object_type => 'ProteinAlignFeature',
      display     => $display,
      renderers   => $self->{'alignment_renderers'},
      strand      => 'b',
    });
  }
}

sub add_protein_features {
  my ($self, $key, $hashref) = @_;

  # We have two separate glyphsets in this in this case
  # P_feature and P_domain - plus domains get copied onto gsv_domain as well
  my %menus = (
    domain     => [ 'domain',    'P_domain',   'normal' ],
    feature    => [ 'feature',   'P_feature',  'normal' ],
    alignment  => [ 'alignment', 'P_domain',   'off'    ],
    gsv_domain => [ 'domain',    'gsv_domain', 'normal' ]
  );

  return unless grep $self->get_node($_), keys %menus;

  my ($keys, $data) = $self->_merge($hashref->{'protein_feature'});

  foreach my $menu_code (keys %menus) {
    my $menu = $self->get_node($menu_code);
    
    next unless $menu;
    
    my $type     = $menus{$menu_code}[0];
    my $gset     = $menus{$menu_code}[1];
    my $renderer = $menus{$menu_code}[2];
    
    foreach (@$keys) {
      next if $self->tree->get_node("${type}_$_");
      next if $type ne ($data->{$_}{'type'} || 'feature'); # Don't separate by db in this case
      
      $self->generic_add($menu, $key, "${type}_$_", $data->{$_}, {
        glyphset  => $gset,
        colourset => 'protein_feature',
        display   => $renderer,
        depth     => 1e6,
        strand    => $gset =~ /P_/ ? 'f' : 'b',
      });
    }
  }
}

sub add_repeat_features {
  my ($self, $key, $hashref) = @_;
  my $menu = $self->get_node('repeat');
  
  return unless $menu && $hashref->{'repeat_feature'}{'rows'} > 0;
  
  my $data    = $hashref->{'repeat_feature'}{'analyses'};
  my %options = (
    glyphset    => '_repeat',
    optimizable => 1,
    depth       => 0.5,
    bump_width  => 0,
    strand      => 'r',
  );
  
  $menu->append($self->create_track("repeat_$key", 'All repeats', {
    db          => $key,
    logic_names => [ undef ], # All logic names
    types       => [ undef ], # All repeat types
    name        => 'All repeats',
    description => 'All repeats',
    colourset   => 'repeat',
    display     => 'off',
    renderers   => [qw(off Off normal On)],
    %options
  }));
  
  my $flag    = keys %$data > 1;
  my $colours = $self->species_defs->colour('repeat');
  
  foreach my $key_2 (sort { $data->{$a}{'name'} cmp $data->{$b}{'name'} } keys %$data) {
    if ($flag) {
      # Add track for each analysis
      $self->generic_add($menu, $key, "repeat_${key}_$key_2", $data->{$key_2}, {
        logic_names => [ $key_2 ], # Restrict to a single supset of logic names
        types       => [ undef  ],
        colours     => $colours,
        description => $data->{$key_2}{'desc'},
        display     => 'off',
        %options
      });
    }
    
    my $d2 = $data->{$key_2}{'types'};
    
    if (keys %$d2 > 1) {
      foreach my $key_3 (sort keys %$d2) {
        my $n  = $key_3;
           $n .= " ($data->{$key_2}{'name'})" unless $data->{$key_2}{'name'} eq 'Repeats';
         
        # Add track for each repeat_type;        
        $menu->append($self->create_track('repeat_' . $key . '_' . $key_2 . '_' . $key_3, $n, {
          db          => $key,
          logic_names => [ $key_2 ],
          types       => [ $key_3 ],
          name        => $n,
          colours     => $colours,
          description => "$data->{$key_2}{'desc'} ($key_3)",
          display     => 'off',
          renderers   => [qw(off Off normal On)],
          %options
        }));
      }
    }
  }
}

sub add_simple_features {
  my ($self, $key, $hashref) = @_;
  my $menu = $self->get_node('simple');
  
  return unless $menu;
  
  my ($keys, $data) = $self->_merge($hashref->{'simple_feature'});
  
  foreach (grep !$data->{$_}{'transcript_associated'}, @$keys) {  
    # Allow override of default glyphset, menu etc.
    $menu = $self->get_node($data->{$_}{'menu'}) if $data->{$_}{'menu'};
    
    next unless $menu;
    
    my $glyphset = $data->{$_}{'glyphset'} ? $data->{$_}{'glyphset'}: '_simple';
    my %options  = (
      glyphset  => $glyphset,
      colourset => 'simple',
      strand    => 'r',
    );

    foreach my $opt ('renderers', 'height') {
      $options{$opt} = $data->{$_}{$opt} if $data->{$_}{$opt};
    }
    
    $self->generic_add($menu, $key, "simple_${key}_$_", $data->{$_}, \%options);
  }
}

sub add_decorations {
  my ($self, $key, $hashref) = @_;
  my $menu = $self->get_node('decorations');
  
  return unless $menu;
  
  if ($key eq 'core' && $hashref->{'karyotype'}{'rows'} > 0 && !$self->get_node('ideogram')) {
    $menu->append($self->create_track("chr_band_$key", 'Chromosome bands', {
      db          => $key,
      glyphset    => 'chr_band',
      display     => 'normal',
      strand      => 'f',
      description => 'Cytogenetic bands',
      colourset   => 'ideogram',
      sortable    => 1,
    }));
  }
  
  if ($key eq 'core' && $hashref->{'assembly_exception'}{'rows'} > 0) {
    $menu->append($self->create_track("assembly_exception_$key", 'Assembly exceptions', {
      db           => $key,
      glyphset     => 'assemblyexception',
      height       => 2,
      display      => 'collapsed',
      renderers    => [ 'off', 'Off', 'collapsed', 'Collapsed', 'collapsed_label', 'Collapsed with labels', 'normal', 'Expanded' ],
      strand       => 'x',
      label_strand => 'r',
      short_labels => 0,
      description  => 'GRC assembly patches, haplotype (HAPs) and pseudo autosomal regions (PARs)',
      colourset    => 'assembly_exception',
    }));
  }
  
  if ($key eq 'core' && $hashref->{'misc_feature'}{'sets'}{'NoAnnotation'}) {
    $menu->append($self->create_track('annotation_status', 'Annotation status', {
      db            => $key,
      glyphset      => 'annotation_status',
      height        => 2,
      display       => 'normal',
      strand        => 'x',
      label_strand  => 'r',
      short_labels  => 0,
      depth         => 0,
      description   => 'Unannotated regions',
      colourset     => 'annotation_status',
    }));
  }
}

#----------------------------------------------------------------------#
# Functions to add tracks from compara like databases                  #
#----------------------------------------------------------------------#

sub add_synteny {
  my ($self, $key, $hashref, $species) = @_;
  my $menu = $self->get_node('synteny');
  
  return unless $menu;
  
  my @synteny_species = sort keys %{$hashref->{'SYNTENY'}{$species} || {}};

  return unless @synteny_species;
  
  my $species_defs = $self->species_defs;
  my $colours      = $species_defs->colour('synteny');
  my $self_label   = $species_defs->species_label($species, 'no_formatting');

  foreach my $species_2 (@synteny_species) {
    (my $species_readable = $species_2) =~ s/_/ /g;
    my ($a, $b) = split / /, $species_readable;
    my $caption = substr($a, 0, 1) . ".$b synteny";
    my $label   = $species_defs->species_label($species_2, 'no_formatting');
    (my $name   = "Synteny with $label") =~ s/<.*?>//g;
    
    $menu->append($self->create_track("synteny_$species_2", $name, {
      db          => $key,
      glyphset    => '_synteny',
      species     => $species_2,
      species_hr  => $species_readable,
      caption     => $caption,
      description => qq{<a href="/info/docs/compara/analyses.html#synteny" class="cp-external">Synteny regions</a> between $self_label and $label},
      colours     => $colours,
      display     => 'off',
      renderers   => [qw(off Off normal On)],
      height      => 4,
      strand      => 'r',
    }));
  }
}

sub add_alignments {
  my ($self, $key, $hashref, $species) = @_;
  
  return unless grep $self->get_node($_), qw(multiple_align pairwise_tblat pairwise_blastz pairwise_other conservation);
  
  my $species_defs = $self->species_defs;
  
  return if $species_defs->ENSEMBL_SITETYPE eq 'Pre';
  
  my $alignments = {};
  my $self_label = $species_defs->species_label($species, 'no_formatting');
  my $static     = $species_defs->ENSEMBL_SITETYPE eq 'Vega' ? '/info/data/comparative_analysis.html' : '/info/docs/compara/analyses.html';
 
  foreach my $row (values %{$hashref->{'ALIGNMENTS'}}) {
    next unless $row->{'species'}{$species};
    
    if ($row->{'class'} =~ /pairwise_alignment/) {
      my ($other_species) = grep { !/^$species$|ancestral_sequences$/ } keys %{$row->{'species'}};
         $other_species ||= $species if scalar keys %{$row->{'species'}} == 1;
      my $other_label     = $species_defs->species_label($other_species, 'no_formatting');
      my ($menu_key, $description, $type);
      
      if ($row->{'type'} =~ /(B?)LASTZ_(\w+)/) {
        next if $2 eq 'PATCH';
        
        $menu_key    = 'pairwise_blastz';
        $type        = sprintf '%sLASTz %s', $1, lc $2;
        $description = "$type pairwise alignments";
      } elsif ($row->{'type'} =~ /TRANSLATED_BLAT/) {
        $type        = '';
        $menu_key    = 'pairwise_tblat';
        $description = 'Trans. BLAT net pairwise alignments';
      } else {
        $type        = ucfirst lc $row->{'type'};
        $type        =~ s/\W/ /g;
        $menu_key    = 'pairwise_align';
        $description = 'Pairwise alignments';
      }
      
      $description  = qq{<a href="$static" class="cp-external">$description</a> between $self_label and $other_label"};
      $description .= " $1" if $row->{'name'} =~ /\((on.+)\)/;

      $alignments->{$menu_key}{$row->{'id'}} = {
        db                         => $key,
        glyphset                   => '_alignment_pairwise',
        name                       => $other_label . ($type ?  " - $type" : ''),
        caption                    => $other_label,
        type                       => $row->{'type'},
        species                    => $other_species,
        method_link_species_set_id => $row->{'id'},
        description                => $description,
        order                      => $other_label,
        colourset                  => 'pairwise',
        strand                     => 'r',
        display                    => 'off',
        renderers                  => [ 'off', 'Off', 'compact', 'Compact', 'normal', 'Normal' ],
      };
    } else {
      my $n_species = grep { $_ ne 'ancestral_sequences' } keys %{$row->{'species'}};
      
      my %options = (
        db                         => $key,
        glyphset                   => '_alignment_multiple',
        short_name                 => $row->{'name'},
        type                       => $row->{'type'},
        species_set_id             => $row->{'species_set_id'},
        method_link_species_set_id => $row->{'id'},
        class                      => $row->{'class'},
        colourset                  => 'multiple',
        strand                     => 'f',
      );
      
      if ($row->{'conservation_score'}) {
        my ($program) = $hashref->{'CONSERVATION_SCORES'}{$row->{'conservation_score'}}{'type'} =~ /(.+)_CONSERVATION_SCORE/;
        
        $options{'description'} = qq{<a href="/info/docs/compara/analyses.html#conservation">$program conservation scores</a> based on the $row->{'name'}};
        
        $alignments->{'conservation'}{"$row->{'id'}_scores"} = {
          %options,
          conservation_score => $row->{'conservation_score'},
          name               => "Conservation score for $row->{'name'}",
          caption            => "$n_species way $program scores",
          order              => sprintf('%12d::%s::%s', 1e12-$n_species*10, $row->{'type'}, $row->{'name'}),
          display            => 'off',
          renderers          => [ 'off', 'Off', 'tiling', 'Tiling array' ],
        };
        
        $alignments->{'conservation'}{"$row->{'id'}_constrained"} = {
          %options,
          constrained_element => $row->{'constrained_element'},
          name                => "Constrained elements for $row->{'name'}",
          caption             => "$n_species way $program elements",
          order               => sprintf('%12d::%s::%s', 1e12-$n_species*10+1, $row->{'type'}, $row->{'name'}),
          display             => 'off',
          renderers           => [ 'off', 'Off', 'compact', 'On' ],
        };
      }
      
      $alignments->{'multiple_align'}{$row->{'id'}} = {
        %options,
        name        => $row->{'name'},
        caption     => $row->{'name'},
        order       => sprintf('%12d::%s::%s', 1e12-$n_species*10-1, $row->{'type'}, $row->{'name'}),
        display     => 'off',
        renderers   => [ 'off', 'Off', 'compact', 'On' ],
        description => qq{<a href="/info/docs/compara/analyses.html#conservation">$n_species way whole-genome multiple alignments</a>.; } . 
                       join('; ', sort map { $species_defs->species_label($_, 'no_formatting') } grep { $_ ne 'ancestral_sequences' } keys %{$row->{'species'}}),
      };
    } 
  }
  
  foreach my $menu_key (keys %$alignments) {
    my $menu = $self->get_node($menu_key);
    next unless $menu;
    
    foreach my $key_2 (sort { $alignments->{$menu_key}{$a}{'order'} cmp  $alignments->{$menu_key}{$b}{'order'} } keys %{$alignments->{$menu_key}}) {
      my $row = $alignments->{$menu_key}{$key_2};
      $menu->append($self->create_track("alignment_${key}_$key_2", $row->{'caption'}, $row));
    }
  }
}

#----------------------------------------------------------------------#
# Functions to add tracks from functional genomics like database       #
#----------------------------------------------------------------------#

# needs configuring so tracks only display if data in species fg_database
sub add_regulation_features {
  my ($self, $key, $hashref) = @_;  
  my $menu = $self->get_node('functional');
  
  return unless $menu;
  
  my $reg_regions       = $menu->append($self->create_submenu('functional_other_regulatory_regions', 'Other regulatory regions'));
  my $methylation_menu  = $reg_regions->before($self->create_submenu('functional_dna_methylation', 'DNA Methylation'));
  my ($keys_1, $data_1) = $self->_merge($hashref->{'feature_set'});
  my ($keys_2, $data_2) = $self->_merge($hashref->{'result_set'});
  my %fg_data           = (%$data_1, %$data_2);
  
  foreach my $key_2 (sort grep { !/Regulatory_Build|seg_/ } @$keys_1, @$keys_2) {
    my $type = $fg_data{$key_2}{'type'};
    
    next if !$type || $type eq 'ctcf';
    
    my @renderers;
    
    if ($fg_data{$key_2}{'renderers'}) {
      push @renderers, $_, $fg_data{$key_2}{'renderers'}{$_} for sort keys %{$fg_data{$key_2}{'renderers'}}; 
    } else {
      @renderers = qw(off Off normal On);
    }
    
    $reg_regions->append($self->create_track("${type}_${key}_$key_2", $fg_data{$key_2}{'name'}, { 
      db          => $key,
      glyphset    => $type,
      sources     => 'undef',
      strand      => 'r',
      labels      => 'on',
      depth       => $fg_data{$key_2}{'depth'}     || 0.5,
      colourset   => $fg_data{$key_2}{'colourset'} || $type,
      display     => $fg_data{$key_2}{'display'}   || 'off', 
      description => $fg_data{$key_2}{'description'},
      logic_name  => $fg_data{$key_2}{'logic_names'}[0],
      renderers   => \@renderers, 
    }));
    
    if ($fg_data{$key_2}{'description'} =~ /cisRED/) {
      $reg_regions->append($self->create_track("${type}_${key}_search", 'cisRED Search Regions', {
        db          => $key,
        glyphset    => 'regulatory_search_regions',
        sources     => 'undef',
        strand      => 'r',
        labels      => 'on',
        depth       => 0.5,
        colourset   => 'regulatory_search_regions',
        description => 'cisRED Search Regions',
        display     => 'off',
      }));
    }
  }
  
  # Add internal methylation tracks
  my $db_tables   = $self->databases->{'DATABASE_FUNCGEN'}{'tables'};
  my $methylation = $db_tables->{'methylation'};
  
  foreach my $k (sort { $methylation->{$a}{'description'} cmp $methylation->{$b}{'description'} } keys %$methylation) {
    $methylation_menu->append($self->create_track("methylation_$k", $methylation->{$k}{'name'}, {
      data_id      => $k,
      description  => $methylation->{$k}{'description'},
      strand       => 'r',
      nobump       => 1,
      addhiddenbgd => 1,
      display      => 'off',
      renderers    => [ qw(off Off compact On) ],
      glyphset     => 'fg_methylation',
      colourset    => 'seq',
    }));
  }
  
  $self->add_track('information', 'fg_methylation_legend', 'Methylation Legend', 'fg_methylation_legend', { strand => 'r' });        
}

sub add_regulation_builds {
  my ($self, $key, $hashref) = @_;  
  my $menu = $self->get_node('functional');
  
  return unless $menu;
  
  my ($keys_1, $data_1) = $self->_merge($hashref->{'feature_set'});
  my ($keys_2, $data_2) = $self->_merge($hashref->{'result_set'});
  my %fg_data           = (%$data_1, %$data_2);
  my $key_2             = 'Regulatory_Build';
  my $type              = $fg_data{$key_2}{'type'};
  
  return unless $type;
  
  my $hub = $self->hub;
  my $db  = $hub->database('funcgen', $self->species);
  
  return unless $db;

  $menu = $menu->append($self->create_submenu('regulatory_features', 'Regulatory features'));
  
  my $db_tables     = $self->databases->{'DATABASE_FUNCGEN'}{'tables'};
  my $reg_feats     = $menu->append($self->create_submenu('reg_features', 'Regulatory features'));
  my $reg_segs      = $menu->append($self->create_submenu('seg_features', 'Segmentation features'));
  my $adaptor       = $db->get_FeatureTypeAdaptor;
  my $evidence_info = $adaptor->get_regulatory_evidence_info;
  my @cell_lines    = sort { ($b eq 'MultiCell') <=> ($a eq 'MultiCell') || $a cmp $b } map [split ':']->[0], keys %{$db_tables->{'cell_type'}{'ids'}}; # Put MultiCell first
  my (@renderers, $prev_track, %matrix_menus, %matrix_rows);
  
  # FIXME: put this in db
  my %default_evidence_types = (
    CTCF     => 1,
    DNase1   => 1,
    H3K4me3  => 1,
    H3K36me3 => 1,
    H3K27me3 => 1,
    H3K9me3  => 1,
    PolII    => 1,
    PolIII   => 1,
  );
  
  if ($fg_data{$key_2}{'renderers'}) {
    push @renderers, $_, $fg_data{$key_2}{'renderers'}{$_} for sort keys %{$fg_data{$key_2}{'renderers'}}; 
  } else {
    @renderers = qw(off Off normal On);
  }
  
  foreach my $cell_line (@cell_lines) { 
    ### Add tracks for cell_line peaks and wiggles only if we have data to display
    my $ftypes     = $db_tables->{'regbuild_string'}{'feature_type_ids'}{$cell_line}      || {};
    my $focus_sets = $db_tables->{'regbuild_string'}{'focus_feature_set_ids'}{$cell_line} || {};
    my @sets;
    
    push @sets, 'core'     if scalar keys %$focus_sets && scalar keys %$focus_sets <= scalar keys %$ftypes;
    push @sets, 'non_core' if scalar keys %$ftypes != scalar keys %$focus_sets && $cell_line ne 'MultiCell';
    
    foreach my $set (@sets) {
      $matrix_menus{$set} ||= [ "reg_feats_$set", $evidence_info->{$set}{'name'}, {
        menu   => 'matrix',
        url    => $hub->url('Config', { action => 'Matrix', function => $hub->action, partial => 1, menu => "reg_feats_$set" }),
        matrix => {
          section     => $menu->data->{'caption'},
          header      => $evidence_info->{$set}{'long_name'},
          description => $db_tables->{'feature_set'}{'analyses'}{'Regulatory_Build'}{'desc'}{$set},
          axes        => { x => 'Cell type', y => 'Evidence type' },
        }
      }];
      
      foreach (@{$evidence_info->{$set}{'classes'}}) {
        foreach (@{$adaptor->fetch_all_by_class($_)}) {
          $matrix_rows{$cell_line}{$set}{$_->name} ||= { row => $_->name, group => $_->class, group_order => $_->class =~ /^(Polymerase|Open Chromatin)$/ ? 1 : 2, on => $default_evidence_types{$_->name} } if $ftypes->{$_->dbID};
        }
      }
    }
  }
  
  $matrix_menus{$_} = $menu->after($self->create_submenu(@{$matrix_menus{$_}})) for 'non_core', 'core';
  
  foreach my $cell_line (@cell_lines) {
    my $track_key = "reg_feats_$cell_line";
    my $display   = 'off';
    my ($label, %evidence_tracks);
    
    if ($cell_line eq 'MultiCell') {  
      $display = $fg_data{$key_2}{'display'} || 'off';
    } else {
      $label = ": $cell_line";
    }
    
    $prev_track = $reg_feats->append($self->create_track($track_key, "$fg_data{$key_2}{'name'}$label", {
      db          => $key,
      glyphset    => $type,
      sources     => 'undef',
      strand      => 'r',
      depth       => $fg_data{$key_2}{'depth'}     || 0.5,
      colourset   => $fg_data{$key_2}{'colourset'} || $type,
      description => $fg_data{$key_2}{'description'}{'reg_feats'},
      display     => $display,
      renderers   => \@renderers,
      cell_line   => $cell_line
    }));
    
    if ($fg_data{"seg_$cell_line"}{'key'} eq "seg_$cell_line") {
      $prev_track = $reg_segs->append($self->create_track("seg_$cell_line", "Reg. Segs: $cell_line", {
        db          => $key,
        glyphset    => 'fg_segmentation_features',
        sources     => 'undef',
        strand      => 'r',
        labels      => 'on',
        depth       => 0,
        colourset   => 'fg_segmentation_features',
        display     => 'off',
        description => $fg_data{"seg_$cell_line"}{'description'},
        renderers   => \@renderers,
        cell_line   => $cell_line,
        caption     => "Reg. Segs $cell_line",
        track_after => $prev_track,
      }));
    }
    
    my %column_data = (
      db        => $key,
      glyphset  => 'fg_multi_wiggle',
      strand    => 'r',
      depth     => $fg_data{$key_2}{'depth'} || 0.5,
      colourset => 'feature_set',
      cell_line => $cell_line,
      menu_key  => 'regulatory_features',
      renderers => [
        'off',            'Off', 
        'compact',        'Peaks', 
        'tiling',         'Signal', 
        'tiling_feature', 'Both' 
      ],
    );
    
    foreach (grep exists $matrix_rows{$cell_line}{$_}, keys %matrix_menus) {
      $prev_track = $self->add_matrix({
        track_name  => "$evidence_info->{$_}{'name'}$label",
        track_after => $prev_track,
        matrix      => {
          menu   => $matrix_menus{$_}->id,
          column => $cell_line,
          rows   => [ values %{$matrix_rows{$cell_line}{$_}} ],
        },
        column_data => {
          set         => $_,
          label       => "$evidence_info->{$_}{'label'} $label",
          description => $fg_data{$key_2}{'description'}{$_},
          %column_data
        }, 
      }, $matrix_menus{$_});
    }
  }
  
  if ($db_tables->{'cell_type'}{'ids'}) {
    $self->add_track('information', 'fg_regulatory_features_legend',   'Reg. Features Legend',          'fg_regulatory_features_legend',   { strand => 'r', colourset => 'fg_regulatory_features'   });        
    $self->add_track('information', 'fg_segmentation_features_legend', 'Reg. Segments Legend',          'fg_segmentation_features_legend', { strand => 'r', colourset => 'fg_segmentation_features' });
    $self->add_track('information', 'fg_multi_wiggle_legend',          'Cell/Tissue Regulation Legend', 'fg_multi_wiggle_legend',          { strand => 'r', display => 'off' });
  }
}

sub add_oligo_probes {
  my ($self, $key, $hashref) = @_; 
  my $menu = $self->get_node('oligo');
  
  return unless $menu;
  
  my $data        = $hashref->{'oligo_feature'}{'arrays'};
  my $description = $hashref->{'oligo_feature'}{'analyses'}{'AlignAffy'}{'desc'};  # Different loop - no analyses - base on probeset query results
  
  foreach my $key_2 (sort keys %$data) {
    my $key_3 = $key_2; 
    $key_2    =~ s/:/__/;  
    
    $menu->append($self->create_track("oligo_${key}_" . uc $key_2, $key_3, {
      glyphset    => '_oligo',
      db          => $key,
      sub_type    => 'oligo',
      array       => $key_2,
      object_type => 'ProbeFeature',
      colourset   => 'feature',
      description => $description,
      caption     => $key_3,
      strand      => 'b',
      display     => 'off',
      renderers   => $self->{'alignment_renderers'}
    }));
  }
}

#----------------------------------------------------------------------#
# Functions to add tracks from variation like databases                #
#----------------------------------------------------------------------#

sub add_sequence_variations {
  my ($self, $key, $hashref) = @_;
  my $menu = $self->get_node('variation');
  
  return unless $menu && $hashref->{'variation_feature'}{'rows'} > 0;
  
  my $options = {
    db         => $key,
    glyphset   => '_variation',
    strand     => 'r',
    depth      => 0.5,
    bump_width => 0,
    colourset  => 'variation',
    display    => 'off',
    renderers  => [ 'off', 'Off', 'normal', 'Normal (collapsed for windows over 200kb)', 'compact', 'Collapsed', 'labels', 'Expanded with name (hidden for windows over 10kb)', 'nolabels', 'Expanded without name' ],
  };
  
  if (defined($hashref->{'menu'}) && scalar @{$hashref->{'menu'}} && grep {$_->{key} =~ /dbsnp/i} @{$hashref->{menu}}) {
    $self->add_sequence_variations_meta($key, $hashref, $options);
  } else {
    $self->add_sequence_variations_default($key, $hashref, $options);
  }

  $self->add_track('information', 'variation_legend', 'Variation Legend', 'variation_legend', { strand => 'r' });
}

# adds variation tracks in structure defined in variation meta table
sub add_sequence_variations_meta {
  my ($self, $key, $hashref, $options) = @_;
  my $menu = $self->get_node('variation');
  
  foreach my $menu_item(sort {$a->{type} cmp $b->{type} || $a->{parent} cmp $b->{parent}} @{$hashref->{'menu'}}) {
    next if $menu_item->{'type'} =~  /^sv_/; # exclude structural variant items
    
    my $node;
    
    if ($menu_item->{'type'} eq 'menu') { # just a named submenu
      $node = $self->create_submenu($menu_item->{'key'}, $menu_item->{'long_name'});
    } elsif ($menu_item->{'type'} eq 'source') { # source type
      (my $temp_name     = $menu_item->{'long_name'}) =~ s/ variants$//;
      my $other_sources = ($menu_item->{'long_name'} =~ /all other sources/);
      
      $node = $self->create_track($menu_item->{'key'}, $menu_item->{'long_name'}, {
        %$options,
        caption     => $menu_item->{'long_name'},
        sources     => $other_sources ? undef : [ $temp_name ],
        description => $other_sources ? 'Sequence variants from all sources' : $hashref->{'source'}{'descriptions'}{$temp_name},
      });
    } elsif ($menu_item->{'type'} eq 'set') { # set type
      (my $temp_name = $menu_item->{'key'})       =~ s/^variation_set_//;
      (my $caption   = $menu_item->{'long_name'}) =~ s/1000 Genomes/1KG/;  # shorten name for side of image
      (my $set_name  = $menu_item->{'long_name'}) =~ s/All HapMap/HapMap/; # hack for HapMap set name - remove once variation team fix data for 68
      
      next if $set_name =~ /HapMap.+/;
      
      $node = $self->create_track($menu_item->{'key'}, $menu_item->{'long_name'}, {
        %$options,
        caption     => $caption,
        sources     => undef,
        sets        => [ $set_name ],
        set_name    => $set_name,
        description => $hashref->{'variation_set'}{'descriptions'}{$temp_name}
      });
    }
    
    # get the node onto which we're going to add this item, then append it
    ($self->get_node($menu_item->{'parent'}) || $menu)->append($node) if $node;
  }
}

# adds variation tracks the old, hacky way
sub add_sequence_variations_default {
  my ($self, $key, $hashref, $options) = @_;
  my $menu               = $self->get_node('variation');
  my $sequence_variation = $self->create_submenu('sequence_variations', 'Sequence variants');

  # XXX: This is a hack until Laurent gets meta information into the
  #        private variation database. It should be there for e76, at
  #        which point this code should never by called in that context.
  my $title = 'Sequence variants (all sources)';
  $title = 'Sequence variants (DECIPHER/LOVD)' if $key eq 'variation_private';

  $sequence_variation->append($self->create_track("variation_feature_$key", $title, {
    %$options,
    sources     => undef,
    description => 'Sequence variants from all sources',
  }));
  
  foreach my $key_2 (sort keys %{$hashref->{'source'}{'counts'} || {}}) {
    next unless $hashref->{'source'}{'counts'}{$key_2} > 0;
    next if     $hashref->{'source'}{'somatic'}{$key_2} == 1;
    
    $sequence_variation->append($self->create_track("variation_feature_${key}_$key_2", "$key_2 variations", {
      %$options,
      caption     => $key_2,
      sources     => [ $key_2 ],
      description => $hashref->{'source'}{'descriptions'}{$key_2},
    }));
  }
  
  $menu->append($sequence_variation);

  # add in variation sets
  if ($hashref->{'variation_set'}{'rows'} > 0) {
    my $variation_sets = $self->create_submenu('variation_sets', 'Variation sets');
    
    $menu->append($variation_sets);
    
    foreach my $toplevel_set (
      sort { !!scalar @{$a->{'subsets'}} <=> !!scalar @{$b->{'subsets'}} } 
      sort { $a->{'name'} =~ /^failed/i  <=> $b->{'name'} =~ /^failed/i  } 
      sort { $a->{'name'} cmp $b->{'name'} } 
      values %{$hashref->{'variation_set'}{'supersets'}}
    ) {
      my $name          = $toplevel_set->{'name'};
      my $caption       = $name . (scalar @{$toplevel_set->{'subsets'}} ? ' (all data)' : '');
      my $key           = $toplevel_set->{'short_name'};
      my $set_variation = scalar @{$toplevel_set->{'subsets'}} ? $self->create_submenu("set_variation_$key", $name) : $variation_sets;
      
      $set_variation->append($self->create_track("variation_set_$key", $caption, {
        %$options,
        caption     => $caption,
        sources     => undef,
        sets        => [ $name ],
        set_name    => $name,
        description => $toplevel_set->{'description'},
      }));
      
      # add in sub sets
      if (scalar @{$toplevel_set->{'subsets'}}) {
        foreach my $subset_id (sort @{$toplevel_set->{'subsets'}}) {
          my $sub_set             = $hashref->{'variation_set'}{'subsets'}{$subset_id};
          my $sub_set_name        = $sub_set->{'name'}; 
          my $sub_set_description = $sub_set->{'description'};
          my $sub_set_key         = $sub_set->{'short_name'};
          
          $set_variation->append($self->create_track("variation_set_$sub_set_key", $sub_set_name, {
            %$options,
            caption     => $sub_set_name,
            sources     => undef,
            sets        => [ $sub_set_name ],
            set_name    => $sub_set_name,
            description => $sub_set_description
          }));
        }
       
        $variation_sets->append($set_variation);
      }
    }
  }
}

sub add_phenotypes {
  my ($self, $key, $hashref) = @_;
  
  return unless $hashref->{'phenotypes'}{'rows'} > 0;
  
  my $p_menu = $self->get_node('phenotype');
  
  unless($p_menu) {
    my $menu = $self->get_node('variation');
    return unless $menu;
    $p_menu = $self->create_submenu('phenotype', 'Phenotype annotations');
    $menu->append($p_menu);
  }
  
  return unless $p_menu;
  
  my $pf_menu = $self->create_submenu('phenotype_features', 'Phenotype annotations');
  
  my %options = (
    db => $key,
    glyphset => 'phenotype_feature',
    depth      => '5',
    bump_width => 0,
    colourset  => 'phenotype_feature',
    display    => 'off',
    strand     => 'r',
    renderers  => [ 'off', 'Off', 'gene_nolabel', 'Expanded', 'compact', 'Compact' ],
  );
  
  $pf_menu->append($self->create_track('phenotype_all', 'Phenotype annotations (all types)', {
    %options,
    caption => 'Phenotypes',
    type => undef,
    description => 'Phenotype annotations on '.(join ", ", map {$_.'s'} keys %{$hashref->{'phenotypes'}{'types'}}),
  }));
  
  foreach my $type(keys %{$hashref->{'phenotypes'}{'types'}}) {
    $pf_menu->append($self->create_track('phenotype_'.lc($type), 'Phenotype annotations ('.$type.'s)', {
      %options,
      caption => 'Phenotypes ('.$type.'s)',
      type => $type,
      description => 'Phenotype annotations on '.$type.'s',
    }));
  }
  
  $p_menu->append($pf_menu);
}

sub add_structural_variations {
  my ($self, $key, $hashref) = @_;
  my $menu = $self->get_node('variation');
  my @A = keys $hashref;
  
  return unless $menu && scalar(keys(%{$hashref->{'structural_variation'}{'counts'}})) > 0;
  
  my $sv_menu             = $self->create_submenu('structural_variation', 'Structural variants');
  my $structural_variants = $self->create_submenu('structural_variants',  'Structural variants');
  my $desc                = 'The colours correspond to the structural variant classes.';
     $desc               .= '<br />For an explanation of the display, see the <a rel="external" href="http://www.ncbi.nlm.nih.gov/dbvar/content/overview/#representation">dbVar documentation</a>.';
  my %options             = (
    glyphset   => 'structural_variation',
    strand     => 'r', 
    bump_width => 0,
    height     => 6,
    depth      => 100,
    colourset  => 'structural_variant',
    display    => 'off',
    renderers  => [ 'off', 'Off', 'compact', 'Compact', 'gene_nolabel', 'Expanded' ],
  );
  
  # Complete overlap (Larger structural variants)
  $structural_variants->append($self->create_track('variation_feature_structural_larger', 'Larger structural variants (all sources)', {   
    %options,
    db         => 'variation',
    caption     => 'Larger structural variants',
    sources     => undef,
    description => "Structural variants from all sources which are at least 1Mb in length. $desc",
    min_size    => 1e6,
  }));
  
  # Partial overlap (Smaller structural variants)
  $structural_variants->append($self->create_track('variation_feature_structural_smaller', 'Smaller structural variants (all sources)', {   
    %options,
    db         => 'variation',
    caption     => 'Smaller structural variants',
    sources     => undef,
    description => "Structural variants from all sources which are less than 1Mb in length. $desc",
    depth       => 10,
    max_size    => 1e6 - 1,
  }));
  
  foreach my $key_2 (sort keys %{$hashref->{'structural_variation'}{'counts'} || {}}) {    
    ## FIXME - Nasty hack to get variation tracks correctly configured
    my $db = $key_2 =~ /DECIPHER/ ? 'variation_private' : 'variation';
    $structural_variants->append($self->create_track("variation_feature_structural_$key_2", "$key_2 structural variations", {
      %options,
      db          => $db,
      caption     => $key_2,
      source      => $key_2,
      description => $hashref->{'source'}{'descriptions'}{$key_2},
    }));  
  }
  
  # Structural variation sets and studies
  foreach my $menu_item (@{$hashref->{'menu'}}) {
    next if $menu_item->{'type'} !~ /^sv_/;
    
    my $node_name = "$menu_item->{'long_name'} (structural variants)";
    
    if ($menu_item->{'type'} eq 'sv_set') {
      my $temp_name = $menu_item->{'key'};
         $temp_name =~ s/^sv_set_//;
      
      $structural_variants->append($self->create_track($menu_item->{'key'}, $node_name, {
        %options,
        caption     => $node_name,
        sources     => undef,
        sets        => [ $menu_item->{'long_name'} ],
        set_name    => $menu_item->{'long_name'},
        description => $hashref->{'variation_set'}{'descriptions'}{$temp_name},
      }));
    }
    elsif ($menu_item->{'type'} eq 'sv_study') {
      my $name = $menu_item->{'key'};
      
      $structural_variants->append($self->create_track($name, $node_name, {
        %options,
        caption     => $node_name,
        sources     => undef,
        study       => [ $name ],
        study_name  => $name,
        description => 'DGVa study: '.$hashref->{'structural_variation'}{'study'}{'descriptions'}{$name},
      }));  
    }
  }
  
  $self->add_track('information', 'structural_variation_legend', 'Structural Variation Legend', 'structural_variation_legend', { strand => 'r' });
  
  $sv_menu->append($structural_variants);
  $menu->append($sv_menu);
}
  
sub add_copy_number_variant_probes {
  my ($self, $key, $hashref) = @_;
  my $menu = $self->get_node('variation');
  
  return unless $menu && scalar(keys(%{$hashref->{'structural_variation'}{'cnv_probes'}{'counts'}})) > 0;
  
  my $sv_menu        = $self->get_node('structural_variation') || $menu->append($self->create_submenu('structural_variation', 'Structural variants'));
  my $cnv_probe_menu = $self->create_submenu('cnv_probe','Copy number variant probes');
  
  my %options = (
    db         => $key,
    glyphset   => 'cnv_probes',
    strand     => 'r', 
    bump_width => 0,
    height     => 6,
    colourset  => 'structural_variant',
    display    => 'off'
  );
  
  $cnv_probe_menu->append($self->create_track('variation_feature_cnv', 'Copy number variant probes (all sources)', {   
    %options,
    caption     => 'CNV probes',
    sources     => undef,
    depth       => 10,
    description => 'Copy number variant probes from all sources'
  }));
  
  foreach my $key_2 (sort keys %{$hashref->{'structural_variation'}{'cnv_probes'}{'counts'} || {}}) {   
    $cnv_probe_menu->append($self->create_track("variation_feature_cnv_$key_2", "$key_2", {
      %options,
      caption     => $key_2,
      source      => $key_2,
      depth       => 0.5,
      description => $hashref->{'source'}{'descriptions'}{$key_2}
    }));  
  }
  
  $sv_menu->append($cnv_probe_menu);
}

# The recombination menu contains tracks with information pertaining to variation project, but these tracks actually simple_features stored in the core database
# As core databases are loaded before variation databases, the recombination submenu appears at the top of the variation menu tree, which isn't desirable.
# This function moves it to the end of the tree.
sub add_recombination {
  my ($self, @args) = @_;
  my $menu   = $self->get_node('recombination');
  my $parent = $self->get_node('variation');
  
  $parent->append($menu) if $menu && $parent;
}

sub add_somatic_mutations {
  my ($self, $key, $hashref) = @_;
  
  # check we have any sources with somatic data
  return unless $hashref->{'source'}{'somatic'} && grep {$_} values %{$hashref->{'source'}{'somatic'}};
  
  my $menu = $self->get_node('somatic');
  return unless $menu;
  
  my $somatic = $self->create_submenu('somatic_mutation', 'Somatic variants');
  my %options = (
    db         => $key,
    glyphset   => '_variation',
    strand     => 'r',
    depth      => 0.5,
    bump_width => 0,
    colourset  => 'variation',
    display    => 'off',
    renderers  => [ 'off', 'Off', 'normal', 'Normal (collapsed for windows over 200kb)', 'compact', 'Collapsed', 'labels', 'Expanded with name (hidden for windows over 10kb)', 'nolabels', 'Expanded without name' ],
  );
  
  # All sources
  $somatic->append($self->create_track("somatic_mutation_all", "Somatic variants (all sources)", {
    %options,
    caption     => 'Somatic variants (all sources)',
    description => 'Somatic variants from all sources'
  }));
  
   
  # Mixed source(s)
  foreach my $key_1 (keys(%{$self->species_defs->databases->{'DATABASE_VARIATION'}{'SOMATIC_MUTATIONS'}})) {
    if ($self->species_defs->databases->{'DATABASE_VARIATION'}{'SOMATIC_MUTATIONS'}{$key_1}{'none'}) {
      (my $k = $key_1) =~ s/\W/_/g;
      $somatic->append($self->create_track("somatic_mutation_$k", "$key_1 somatic variants", {
        %options,
        caption     => "$key_1 somatic variants",
        source      => $key_1,
        description => "Somatic variants from $key_1"
      }));
    }
  }
  
  # Somatic source(s)
  foreach my $key_2 (sort grep { $hashref->{'source'}{'somatic'}{$_} == 1 } keys %{$hashref->{'source'}{'somatic'}}) {
    next unless $hashref->{'source'}{'counts'}{$key_2} > 0;
    
    $somatic->append($self->create_track("somatic_mutation_$key_2", "$key_2 somatic mutations (all)", {
      %options,
      caption     => "$key_2 somatic mutations (all)",
      source      => $key_2,
      description => "All somatic variants from $key_2"
    }));
    
    my $tissue_menu = $self->create_submenu('somatic_mutation_by_tissue', 'Somatic variants by tissue');
    
    ## Add tracks for each tumour site
    my %tumour_sites = %{$self->species_defs->databases->{'DATABASE_VARIATION'}{'SOMATIC_MUTATIONS'}{$key_2} || {}};
    
    foreach my $description (sort  keys %tumour_sites) {
      next if $description eq 'none';
      
      my $phenotype_id           = $tumour_sites{$description};
      my ($source, $type, $site) = split /\:/, $description;
      my $formatted_site         = $site;
      $site                      =~ s/\W/_/g;
      $formatted_site            =~ s/\_/ /g;
      
      $tissue_menu->append($self->create_track("somatic_mutation_${key_2}_$site", "$key_2 somatic mutations in $formatted_site", {
        %options,
        caption     => "$key_2 $formatted_site tumours",
        filter      => $phenotype_id,
        description => $description
      }));    
    }
    
    $somatic->append($tissue_menu);
  }
  
  $menu->append($somatic);
}

sub add_somatic_structural_variations {
  my ($self, $key, $hashref) = @_;
  my $menu = $self->get_node('somatic');
  
  return unless $menu && scalar(keys(%{$hashref->{'structural_variation'}{'somatic'}{'counts'}})) > 0;
  
  my $somatic = $self->create_submenu('somatic_structural_variation', 'Somatic structural variants');
  
  my %options = (
    db         => $key,
    glyphset   => 'structural_variation',
    strand     => 'r', 
    bump_width => 0,
    height     => 6,
    colourset  => 'structural_variant',
    display    => 'off',
    renderers  => [ 'off', 'Off', 'compact', 'Compact', 'gene_nolabel', 'Expanded' ],
  );
  
  $somatic->append($self->create_track('somatic_sv_feature', 'Somatic structural variants (all sources)', {   
    %options,
    caption     => 'Somatic structural variants',
    sources     => undef,
    description => 'Somatic structural variants from all sources. For an explanation of the display, see the <a rel="external" href="http://www.ncbi.nlm.nih.gov/dbvar/content/overview/#representation">dbVar documentation</a>. In addition, we display the breakpoints in yellow.',
    depth       => 10
  }));
  
  foreach my $key_2 (sort keys %{$hashref->{'structural_variation'}{'somatic'}{'counts'} || {}}) {
    $somatic->append($self->create_track("somatic_sv_feature_$key_2", "$key_2 somatic structural variations", {
      %options,
      caption     => "$key_2 somatic",
      source      => $key_2,
      description => $hashref->{'source'}{'descriptions'}{$key_2},
      depth       => 100
    }));  
  }
  
  $menu->append($somatic);
}

sub share {
  # Remove anything from user settings that is:
  #   Custom data that the user isn't sharing
  #   A track from a datahub that the user isn't sharing
  #   Not for the species in the image
  # Reduced track order of explicitly ordered tracks if they are after custom tracks which aren't shared
  
  my ($self, %shared_custom_tracks) = @_;
  my $user_settings     = EnsEMBL::Web::Root->deepcopy($self->get_user_settings);
  my $species           = $self->species;
  my $user_data         = $self->get_node('user_data');
  my @unshared_datahubs = grep $_->get('datahub_menu') && !$shared_custom_tracks{$_->id}, @{$self->tree->child_nodes};
  my @user_tracks       = map { $_ ? $_->nodes : () } $user_data;
  my %user_track_ids    = map { $_->id => 1 } @user_tracks;
  my %datahub_tracks    = map { $_->id => [ map $_->id, $_->nodes ] } @unshared_datahubs;
  my %to_delete;
  
  foreach (keys %$user_settings) {
    next if $_ eq 'track_order';
    next if $shared_custom_tracks{$_};
    
    my $node = $self->get_node($_);
    
    $to_delete{$_} = 1 unless $node && $node->parent_node; # delete anything that isn't for this species
    $to_delete{$_} = 1 if $user_track_ids{$_};             # delete anything that isn't shared
  }
  
  foreach (@unshared_datahubs) {
    $to_delete{$_} = 1 for grep $user_settings->{$_}, @{$datahub_tracks{$_->id} || []};  # delete anything for tracks in datahubs that aren't shared
  }
  
  # Reduce track orders if custom tracks aren't shared
  if (scalar keys %to_delete) {
    my %track_ids_to_delete  = map { $_ => 1 } keys %to_delete, map { @{$datahub_tracks{$_->id} || []} } @unshared_datahubs;
    my @removed_track_orders = map { $track_ids_to_delete{$_->id} && $_->{'data'}{'node_type'} eq 'track' ? $_->{'data'}{'order'} : () } @{$self->glyphset_configs};
    
    foreach my $order (values %{$user_settings->{'track_order'}{$species}}) {
      my $i = 0;
      
      for (@removed_track_orders) {
        last if $_ > $order;
        $i++;
      }
      
      $i-- if $i && $removed_track_orders[$i] > $order;
      $order -= $i;
    }
  }
  
  foreach (keys %to_delete) {
    delete $user_settings->{$_};
    delete $user_settings->{'track_order'}{$species}{$_} for $_, "$_.f", "$_.r";
  }
  
  delete $user_settings->{'track_order'}{$_} for grep $_ ne $species, keys %{$user_settings->{'track_order'}};
  
  return $user_settings;
}

1;
