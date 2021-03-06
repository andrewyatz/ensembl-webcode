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

package EnsEMBL::Web::ZMenu;

use strict;

use HTML::Entities qw(encode_entities decode_entities);

use base qw(EnsEMBL::Web::Root);

sub new {
  my ($class, $hub, $object) = @_;
  
  my $self = {
    hub      => $hub,
    object   => $object,
    features => [],
  };
  
  bless $self, $class;
  
  $self->new_feature;
  $self->content;
  
  # stored_entries keeps all entries of all plugins in a hash, keyed by order
  
  foreach my $feature (@{$self->{'features'}}) {
    $feature->{'stored_entries'}{$_->{'order'}} = $_ for @{$feature->{'entries'}};
  }
  
  return $self;
}

sub content {}
sub hub     { return $_[0]{'hub'}; }

sub object {
  my $self = shift;
  $self->{'object'} = shift if @_;
  return $self->{'object'};
}

sub caption {
  my $self = shift;
  $self->{'feature'}{'caption'}      = shift if @_;
  $self->{'feature'}{'caption_link'} = shift if @_;
  return $self->{'feature'}{'caption'};
}

sub highlight {
  my $self = shift;
  $self->{'feature'}{'highlight'} = shift if @_;
  return $self->{'feature'}{'highlight'};
}

sub new_feature {
  my $self = shift;
  
  $self->{'feature'} = {
    entries        => [],
    stored_entries => {},
    order          => 1,
    caption        => '',
    highlight      => '',
  };
  
  push @{$self->{'features'}}, $self->{'feature'};
}

sub click_location { return map { $_[0]->hub->param("click_$_") || () } qw(chr start end); }

sub click_data {
  my $self  = shift;
  my $hub   = $self->hub;
  my @click = $self->click_location;
  
  if (scalar @click == 3) {
    my $image_config = $hub->get_imageconfig($hub->param('config'));
    my $node         = $image_config ? $image_config->get_node($hub->param('track')) : undef;
    my $slice        = $node ? $hub->get_adaptor('get_SliceAdaptor')->fetch_by_region('toplevel', @click) : undef;
    
    return { container => $slice, config => $image_config, my_config => $node } if $slice;
  }
  
  return undef;
}

# When adding an entry you can specify ORDER or POSITION.
# ORDER    is used to set the position of all entries. It is auto generated unless specified,
#          and as such should be set explicitly on ALL entries, or NONE. Any other scenario
#          could result in unexpected behaviour. For these cases use POSITION instead.
# POSITION will insert the entry at that position in the menu. Since it increments all subsequent
#          entries' orders, it should only be used to insert at existing positions.
# It is probably best not to use POSITION and ORDER together, just in case something goes wrong.
sub add_entry {
  my ($self, $entry) = @_;
  
  if ($entry->{'position'}) {
    $_->{'order'}++ for grep $_->{'order'} >= $entry->{'position'}, @{$self->{'feature'}{'entries'}}; # increment order for each entry after the given position
    $entry->{'order'} = $entry->{'position'};
    $self->{'feature'}{'order'}++;
  } else {
    $entry->{'order'} ||= $self->{'feature'}{'order'}++;
  }
  
  push @{$self->{'feature'}{'entries'}}, $entry;
}

sub add_subheader {
  my ($self, $label) = @_;
  $self->add_entry({ type => 'subheader', label_html => $label }) if defined $label;
}

# Can be used from plugins to remove entries from previous plugins' menu.
# Requires a list of entry positions to remove
sub remove_entries {
  my $self = shift;
  delete $self->{'feature'}{'stored_entries'}{$_} for @_;
}

# Generic code to grab hold of an existing entry and modify it
sub modify_entry_by {
  my ($self, $key, $entry) = @_;
  
  foreach (@{$self->{'feature'}{'entries'}}) {
    if ($_->{$key} eq $entry->{$key}) {
      $_ = $entry;
      last;
    }
  }
}

# Delete an entry by its value
sub delete_entry_by_value {
  my ($self, $value) = @_;
  
  foreach my $entry (@{$self->{'feature'}{'entries'}}) {
    foreach (keys %$entry) {
      if ($entry->{$_} eq $value) {
        $entry = undef;
        last;
      }
    }
  }
}

# Delete an entry by type
sub delete_entry_by_type {
  my ($self, $type) = @_;
  
  foreach (@{$self->{'feature'}{'entries'}}) {
    if ($_->{'type'} eq $type) {
      $_ = undef;
      last;
    }
  }
}

# Build and print the JSON response
sub render {
  my $self          = shift;
  my %reserved_keys = map { $_ => 1 } qw(type link link_class label label_html order abs_url external update_params); # These keys will not be added to the entries hash as part of %extra
  my @features;
  
  foreach my $feature (@{$self->{'features'}}) {
    my @entries;
    
    foreach (sort { $a <=> $b } keys %{$feature->{'stored_entries'}}) {
      my $entry = $feature->{'stored_entries'}{$_};
      my $type  = encode_entities($entry->{'type'});
      my %extra = map { $reserved_keys{$_} ? () : ($_ => $entry->{$_}) } keys %$entry;
      my $value;
      
      if ($entry->{'link'}) {
        if ($entry->{'abs_url'}) {
          $value = $entry->{'link'};
        } else {
          $value = sprintf(
            '<a href="%s"%s%s>%s%s%s</a>',
            encode_entities(decode_entities($entry->{'link'})), # Decode links before encoding them stops double encoding when the link is created by EnsEMBL::Web::ExtURL->get_url
            $entry->{'external'} ? ' rel="external"' : '',
            $entry->{'link_class'} ? qq{ class="$entry->{'link_class'}"} : '',
            encode_entities($entry->{'label'}),
            $entry->{'label_html'},
            $entry->{'update_params'},
          );
        }
      } else {
        $value = encode_entities($entry->{'label'}) . $entry->{'label_html'};
      }
      
      push @entries, { key => $type, value => $value, %extra };
    }
    
    if (scalar @entries) {
      push @features, {
        caption   => $feature->{'caption_link'} ? $feature->{'caption'} : encode_entities($feature->{'caption'}),
        highlight => $feature->{'highlight'},
        entries   => \@entries,
      };
    }
  }
  
  print $self->jsonify(\@features);
}

1;
