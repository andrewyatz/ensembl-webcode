package EnsEMBL::Web::DOM;

### Serves as a factory for creating Nodes in the DOM tree

use strict;

use base qw(EnsEMBL::Web::Root);

use EnsEMBL::Web::DOM::Node::Element::Generic;

use constant {
  POSSIBLE_HTML_ELEMENTS => [qw(
      a abbr acronym address area b base bdo big blockquote body br button
      caption cite code col colgroup dd del dfn div dl dt em fieldset form
      frame frameset head h1 h2 h3 h4 h5 h6 hr html i iframe img input ins
      kbd label legend li link map meta noframes noscript object ol optgroup
      option p param pre q samp script select small span strong style sub
      sup table tbody td textarea tfoot th thead title tr tt ul var
  )]
};

sub new {
  ## @constructor
  return bless [
    {}, #used classes
    {}  #mapped classes
  ], shift;
}

sub map_element_class {
  ## Maps an element type with a class
  ## When create_element is called with a given node name, the mapped class is instantiated.
  ## @param Hashref of node_name => element_class
  my ($self, $map) = @_;
  $self->[1]{ lc $_ } = $map->{$_} for keys %$map;
}

sub create_document {
  ## Creates document node
  ## @return Node::Document object
  my $self = shift;
  my $node_class = 'EnsEMBL::Web::DOM::Node::Document';
  $self->dynamic_use($node_class);
  return $node_class->new($self);
}

sub create_element {
  ## Creates an element of a given tag name by instantiating the corresponding class
  ## Also adds attributes, inner_HTML/inner_text and child nodes
  ## @param Element type
  ## @param HashRef of {attrib1 => value1, attrib2 => value2} for attributes, inner_HTML/inner_text and children (value for children key should be ref of array as accepted by append_children method).
  ## @return Element subclass object or undef if unknown node_name
  my ($self, $element_name, $attributes)  = @_;

  $element_name = lc $element_name;
  $attributes ||= {};
  
  my $node_class;
  
  # skip 'dynamic_use' of element class if already required once
  unless ($node_class = $self->[0]{$element_name}) {

    $node_class = $self->_get_mapped_element_class($element_name);
    if (!$self->dynamic_use($node_class)) {
      $node_class = undef;
      $_ eq $element_name and $node_class = 'EnsEMBL::Web::DOM::Node::Element::Generic' and last for @{$self->POSSIBLE_HTML_ELEMENTS};
      return unless $node_class;
    }
    $self->[0]{$element_name} = $node_class;
  }

  my $element = $node_class->new($self);
  $element->node_name = $element_name if $node_class eq 'EnsEMBL::Web::DOM::Node::Element::Generic';
  if (exists $attributes->{'inner_HTML'}) {
    $element->inner_HTML(delete $attributes->{'inner_HTML'});
  }
  elsif (exists $attributes->{'inner_text'}) {
    $element->inner_text(delete $attributes->{'inner_text'});
  }
  my @children = @{delete $attributes->{'children'} || []};
  $element->set_attributes($attributes) if scalar keys %$attributes;
  $element->append_children(@children)  if scalar @children;
  return $element;
}

sub create_text_node {
  ## Creates a text node
  ## @param Text string
  ## @return Text node object
  my $self = shift;
  my $node_class = 'EnsEMBL::Web::DOM::Node::Text';
  $self->dynamic_use($node_class);
  my $text_node = $node_class->new($self);
  $text_node->text(shift) if @_;
  return $text_node;
}

sub create_comment {
  ## Creates comment
  ## @param comment string
  ## @return Comment object
  my $self = shift;
  my $node_class = 'EnsEMBL::Web::DOM::Node::Comment';
  $self->dynamic_use($node_class);
  my $comment_node = $node_class->new($self);
  $comment_node->comment(shift) if @_;
  return $comment_node;
}

sub _get_mapped_element_class {
  my ($self, $element_name) = @_;
  return $self->[1]{$element_name} if exists $self->[1]{$element_name};
  return 'EnsEMBL::Web::DOM::Node::Element::Input::'.ucfirst $1 if $element_name =~ /^input([a-z]+)$/;
  return 'EnsEMBL::Web::DOM::Node::Element::'.ucfirst $element_name;
}

1;