# $Id$

package EnsEMBL::Web::Component;

use strict;

use base qw(EnsEMBL::Web::Root Exporter);

use Exporter;

our @EXPORT_OK = qw(cache cache_print);
our @EXPORT = @EXPORT_OK;

use CGI qw(escapeHTML);
use Digest::MD5 qw(md5_hex);
use Text::Wrap qw(wrap);

use Bio::EnsEMBL::DrawableContainer;
use Bio::EnsEMBL::VDrawableContainer;
use Bio::EnsEMBL::ExternalData::DAS::Coordinator;

use EnsEMBL::Web::Document::Image;
use EnsEMBL::Web::Constants;
use EnsEMBL::Web::Form;
use EnsEMBL::Web::RegObj;
use EnsEMBL::Web::TmpFile::Text;

sub new {
  my ($class, $object) = shift;
  
  my $self = {
    'object' => shift,
    'id' => [split /::/, $class]->[-1] . 'Panel'
  };
  
  bless $self, $class;
  $self->_init;
  
  return $self;
}

sub id {
  my $self = shift;
  $self->{'id'} = shift if @_;
  return $self->{'id'};
}

sub object {
  my $self = shift;
  $self->{'object'} = shift if @_;
  return $self->{'object'};
}

sub cacheable {
  my $self = shift;
  $self->{'cacheable'} = shift if @_;
  return $self->{'cacheable'};
}

sub ajaxable {
  my $self = shift;
  $self->{'ajaxable'} = shift if @_;
  return $self->{'ajaxable'};
}

sub configurable {
  my $self = shift;
  $self->{'configurable'} = shift if @_;
  return $self->{'configurable'};
}

sub cache {
  my ($panel, $obj, $type, $name) = @_;
  my $cache = new EnsEMBL::Web::TmpFile::Text(
    prefix   => $type,
    filename => $name,
  );
  return $cache;
}

sub cache_print {
  my ($cache, $string_ref) = @_;
  $cache->print($$string_ref) if $string_ref;
}

sub site_name   { return $SiteDefs::SITE_NAME || $SiteDefs::ENSEMBL_SITETYPE; }
sub image_width { return $ENV{'ENSEMBL_IMAGE_WIDTH'}; }
sub has_image   { return 0; }
sub cache_key   { return undef; }
sub caption     { return undef; }
sub _init       { return; }

sub _error   { return shift->_info_panel('error',   @_);  } # Fatal error message. Couldn't perform action
sub _warning { return shift->_info_panel('warning', @_ ); } # Error message, but not fatal
sub _info    { return shift->_info_panel('info',    @_ ); } # Extra information 
sub _hint    { my ($self, $id, $caption, $desc, $width) = @_; return $self->_info_panel('hint hint_flag', $caption, $desc, $width, $id); } # Extra information, hideable

sub _info_panel {
  my ($self, $class, $caption, $desc, $width, $id) = @_;
  
  return sprintf (
    '<div%s style="width:%s" class="%s"><h3>%s</h3><div class="error-pad">%s</div></div>',
    $id ? qq{ id="$id"} : '',
    $width || $self->image_width . 'px', 
    $class, 
    $caption, 
    $desc
  );
}

sub ajax_url {
  my ($self, $function_name) = @_;
  
  my $object = $self->object;
  my ($ensembl, $plugin, $component, $type, $module) = split '::', ref $self;
  
  my $url = join '/', $object->species_path, 'Component', $object->type, $plugin, $module;
  $url .= "/$function_name" if $function_name && $self->can("content_$function_name");
  $url .= "?$ENV{'QUERY_STRING'}";
  $url .= ';_rmd=' . substr md5_hex($ENV{'REQUEST_URI'}), 0, 4;
  
  return $url;
}

# Attach all das sources from an image config
sub _attach_das {
  my ($self, $wuc) = @_;

  # Look for all das sources which are configured and turned on
  my @das_nodes = map {
    $_->get('glyphset') eq '_das' && $_->get('display') ne 'off' ? @{$_->get('logicnames')||[]} : ()
  }  $wuc->tree->nodes;
  
  return unless @das_nodes; # Return if no sources to be drawn
 
  # Check to see if they really exists, and get entries from get_all_das call
  my %T = %{$ENSEMBL_WEB_REGISTRY->get_all_das($self->object->species)};
  my @das_sources = @T{@das_nodes};
  return unless @das_sources; # Return if no sources exist

  # Cache the DAS Coordinator object (with key das_coord)
  $wuc->cache('das_coord',  
    Bio::EnsEMBL::ExternalData::DAS::Coordinator->new(
      -sources => \@das_sources,
      -proxy   => $self->object->species_defs->ENSEMBL_WWW_PROXY,
      -noproxy => $self->object->species_defs->ENSEMBL_NO_PROXY,
      -timeout => $self->object->species_defs->ENSEMBL_DAS_TIMEOUT
    )
  );
}

# Creates a modal-friendly form with hidden elements to automatically pass 
# _referer and x_requested_with - and to optionally handle wizard buttons
sub modal_form {
  my ($self, $name, $action, $options) = @_;
  
  my $object = $self->object;
  my $form_action = $action;
  
  if ($options->{'wizard'}) {
    my $species = $ENV{'ENSEMBL_TYPE'} eq 'UserData' ? $object->data_species : $object->species;
    
    $form_action  = "/$species" if $species;
    $form_action .= sprintf '/%s/Wizard', $object->type;
  }
  
  my $form = new EnsEMBL::Web::Form($name, $form_action, $options->{'method'} || 'post');
  my $label = $options->{'label'} || 'Next >';

  $form->add_element('type' => 'Hidden', 'name' => '_referer',         'value' => $object->param('_referer'));
  $form->add_element('type' => 'Hidden', 'name' => 'x_requested_with', 'value' => $object->param('x_requested_with'));
  
  if ($options->{'wizard'}) {
    $form->add_button('type' => 'Submit', 'name' => 'wizard_submit', 'value' => '< Back') unless defined $options->{'back_button'} && $options->{'back_button'} == 0;
    
    # Include current and former nodes in _backtrack
    if (my @tracks = $object->param('_backtrack')) {
      foreach my $step (@tracks) {
        next unless $step;
        $form->add_element('type' => 'Hidden', 'name' => '_backtrack', 'value' => $step);
      }
    }
    
    $form->add_button('type'  => 'Submit', 'name' => 'wizard_submit',      'value' => $label);
    $form->add_element('type' => 'Hidden', 'name' => '_backtrack',         'value' => $object->action);
    $form->add_element('type' => 'Hidden', 'name' => 'wizard_next',        'value' => $action);
    $form->add_element('type' => 'Hidden', 'name' => 'wizard_ajax_submit', 'value' => '');
  } elsif (!$options->{'no_button'}) {
    $form->add_button('type' => 'Submit', 'name' => 'submit', 'value' => $label);
  }

  return $form;
}

sub new_image {
  my $self = shift;
  my $object = $self->object;
  
  my %formats = EnsEMBL::Web::Constants::FORMATS;
  my $image_config = ref $_[0] eq 'ARRAY' ? $_[0][1] : $_[1];
  
  $self->id($image_config->{'type'});
  
  # Set text export on image config
  $image_config->set_parameter('text_export', $object->param('export')) if $formats{$object->param('export')}{'extn'} eq 'txt';
  
  my $image = new EnsEMBL::Web::Document::Image($object->species_defs);
  $image->drawable_container = new Bio::EnsEMBL::DrawableContainer(@_);
  $image->prefix($object->prefix) if $object->prefix;
  
  return $image;
}

sub new_vimage {
  my $self = shift;
  my $object = $self->object;
  
  $self->id($_[1]->{'type'}); # $_[1] is image config
  
  my $image = new EnsEMBL::Web::Document::Image($object->species_defs);
  $image->drawable_container = new Bio::EnsEMBL::VDrawableContainer(@_);
  
  return $image;
}

sub new_karyotype_image {
  my ($self, $image_config) = @_;
  my $object = $self->object;
  
  $self->id($image_config->{'type'}) if $image_config;
  
  my $image = new EnsEMBL::Web::Document::Image($object->species_defs);
  $image->{'object'} = $object;
  
  return $image;
}

sub _export_image {
  my ($self, $image, $flag) = @_;
  
  $image->{'export'} = 'iexport' . ($flag ? " $flag" : '');
  
  my ($format, $scale) = $self->object->param('export') ? split /-/, $self->object->param('export'), 2 : ('', 1);
  $scale eq 1 if $scale <= 0;
  
  my %formats = EnsEMBL::Web::Constants::FORMATS;
  
  if ($formats{$format}) {
    $image->drawable_container->{'config'}->set_parameter('sf',$scale);
    (my $comp = ref $self) =~ s/[^\w\.]+/_/g;
    my $filename = sprintf '%s-%s-%s.%s', $comp, $self->object->_filename, $scale, $formats{$format}{'extn'};
    
    if ($self->object->param('download')) {
      $self->object->input->header(-type => $formats{$format}{'mime'}, -attachment => $filename);
    } else {
      $self->object->input->header(-type => $formats{$format}{'mime'}, -inline => $filename);
    }

    if ($formats{$format}{'extn'} eq 'txt') {
      print $image->drawable_container->{'export'};
      return 1;
    }

    $image->render($format);
    return 1;
  }
  
  return 0;
}

sub _matches {
  my ($self, $key, $caption, @keys) = @_;
  
  my $object = $self->object;
  my $label  = $object->species_defs->translate($caption);
  my $obj    = $object->Obj;

  # Check cache
  if (!$object->__data->{'links'}) {
    my @similarity_links = @{$object->get_similarity_hash($obj)};
    
    return unless @similarity_links;
    
    $self->_sort_similarity_links(@similarity_links);
  }

  my @links = map { @{$object->__data->{'links'}{$_}||[]} } @keys;
  
  return unless @links;

  my $db    = $object->get_db;
  my $entry = $object->gene_type || 'Ensembl';

  # add table call here
  my $html;
  
  if ($object->species_defs->ENSEMBL_SITETYPE eq 'Vega') {
    $html = '<p></p>';
  } else {
    $html = "<p><strong>This $entry entry corresponds to the following database identifiers:</strong></p>";
  }
  
  $html .= '<table cellpadding="4">';
  
  @links = $self->remove_redundant_xrefs(@links) if $keys[0] eq 'ALT_TRANS';
  
  return unless @links;
  
  my $old_key = '';
  
  foreach my $link (@links) {
    my ($key, $text) = @$link;
    
    if ($key ne $old_key) {
      $html .= '<div class="small">GO mapping is inherited from swissprot/sptrembl</div>' if $old_key eq 'GO';
      $html .= '</td></tr>' if $old_key ne '';
      $html .= qq{<tr><th style="white-space: nowrap; padding-right: 1em">$key:</th><td>};
      
      $old_key = $key;
    }
    
    $html .= $text;
  }
  
  $html .= '</td></tr></table>';

  return $html;
}

sub _sort_similarity_links {
  my $self = shift;
  my @similarity_links = @_;
  
  my $object   = $self->object;
  my $database = $object->database;
  my $db       = $object->get_db;
  my $urls     = $object->ExtURL;
  my $fv_type  = $object->action eq 'Oligos' ? 'OligoFeature' : 'Xref'; # default link to featureview is to retrieve an Xref
  my (%affy, %exdb);
  
  foreach my $type (sort {
    $b->priority        <=> $a->priority ||
    $a->db_display_name cmp $b->db_display_name ||
    $a->display_id      cmp $b->display_id
  } @similarity_links) {
    my $link = '';
    my $join_links = 0;
    my $externalDB = $type->database;
    my $display_id = $type->display_id;
    my $primary_id = $type->primary_id;
    
    next if $type->status eq 'ORTH';                            # remove all orthologs
    next if lc $externalDB eq 'medline';                        # ditch medline entries - redundant as we also have pubmed
    next if $externalDB =~ /^flybase/i && $display_id =~ /^CG/; # ditch celera genes from FlyBase
    next if $externalDB eq 'Vega_gene';                         # remove internal links to self and transcripts
    next if $externalDB eq 'Vega_transcript';
    next if $externalDB eq 'Vega_translation';
    next if $externalDB eq 'OTTP' && $display_id =~ /^\d+$/;    # don't show vega translation internal IDs
    
    if ($externalDB eq 'GO') {
      push @{$object->__data->{'links'}{'go'}}, $display_id;
      next;
    } elsif ($externalDB eq 'GKB') {
      my ($key, $primary_id) = split ':', $display_id;
      push @{$object->__data->{'links'}{'gkb'}->{$key}}, $type;
      next;
    }
    
    my $text = $display_id;
    
    (my $A = $externalDB) =~ s/_predicted//;
    
    if ($urls and $urls->is_linked($A)) {
      my $link = $urls->get_url($A, $primary_id);
      my $word = $display_id;
      $word .= " ($primary_id)" if $A eq 'MARKERSYMBOL';
      
      if ($link) {
        $text = qq{<a href="$link">$word</a>};
      } else {
        $text = $word;
      }
    }
    
    if ($type->isa('Bio::EnsEMBL::IdentityXref')) {
      $text .= ' <span class="small"> [Target %id: ' . $type->target_identity . '; Query %id: ' . $type->query_identity . ']</span>';
      $join_links = 1;
    }
    
    if ($object->species_defs->ENSEMBL_PFETCH_SERVER && $externalDB =~ /^(SWISS|SPTREMBL|LocusLink|protein_id|RefSeq|EMBL|Gene-name|Uniprot)/i) {
      my $seq_arg = $display_id;
      $seq_arg = "LL_$seq_arg" if $externalDB eq 'LocusLink';
      
      $text .= sprintf ' [<a href="/%s/Transcript/Similarity/Align?t=%s;sequence=%s;db=%s">align</a>] ', $object->species, $object->stable_id, $seq_arg, $db;
    }
    
    $text .= sprintf ' [<a href="%s">Search GO</a>]', $urls->get_url('GOSEARCH', $primary_id) if $externalDB =~ /^(SWISS|SPTREMBL)/i; # add Search GO link;
    
    if ($type->description) {
      (my $D = $type->description) =~ s/^"(.*)"$/$1/;
      
      $text .= '<br />' . escapeHTML($D);
      $join_links = 1;
    }
    
    if ($join_links) {
      $text = qq{\n <div>$text};
    } else {
      $text = qq{\n <div class="multicol">$text};
    }
    
    # override for Affys - we don't want to have to configure each type, and
    # this is an internal link anyway.
    if ($externalDB =~ /^AFFY_/i) {
      next if $affy{$display_id} && $exdb{$type->db_display_name}; # remove duplicates
      
      $text = qq{\n  <div class="multicol"> $display_id};
      $affy{$display_id}++;
      $exdb{$type->db_display_name}++;
    }

    # add link to featureview
    my $link_name = $fv_type eq 'OligoFeature' ? $display_id : $primary_id;
    my $link_type = $fv_type eq 'OligoFeature' ? $fv_type : "${fv_type}_$externalDB";
    
    my $k_url = $object->_url({
      type   => 'Location',
      action => 'Genome',
      id     => $link_name,
      ftype  => $link_type
    });
    
    $text .= qq{  [<a href="$k_url">view all locations</a>]};
    $text .= '</div>';
    
    push @{$object->__data->{'links'}{$type->type}}, [ $type->db_display_name || $externalDB, $text ];
  }
}

sub remove_redundant_xrefs {
  my ($self, @links) = @_;
  my %priorities;

  foreach my $link (@links) {
    my ($key, $text) = @$link;
    $priorities{$key} = $text if $text =~ />OTT|>ENST/;
  }
  
  return () if ref $self->object->Obj eq 'Bio::EnsEMBL::Gene'; # another hack to deal with mouse which has no 'Havana gene'
  
  foreach my $type (
    'Transcript having exact match between ENSEMBL and HAVANA',
    'Ensembl transcript having exact match with Havana',
    'Havana transcript having same CDS',
    'Ensembl transcript sharing CDS with Havana',
    'Havana transcripts'
  ) {
    if ($priorities{$type}) {
      my @munged_links;
      $munged_links[0] = [ $type, $priorities{$type} ];
      return @munged_links;;
    }
  }
  
  return @links;
}

# Simple subroutine to dump a formatted "warn" block to the error logs - useful when debugging complex
# data structures etc... 
# output looks like:
#
#  ###########################
#  #                         #
#  # TEXT. TEXT. TEXT. TEXT. #
#  # TEXT. TEXT. TEXT. TEXT. #
#  # TEXT. TEXT. TEXT. TEXT. #
#  #                         #
#  # TEXT. TEXT. TEXT. TEXT. #
#  # TEXT. TEXT. TEXT. TEXT. #
#  #                         #
#  ###########################
sub _warn_block {
  my $self = shift;
  
  my $width       = 128;
  my $border_char = '#';
  my $template    = sprintf "%s %%-%d.%ds %s\n", $border_char, $width-4,$width-4, $border_char;
  my $line        = $border_char x $width;
  
  warn "\n";
  warn "$line\n";
  
  $Text::Wrap::columns = $width-4;
  
  foreach my $l (@_) {
    my $lines = wrap('','', $l);
    
    warn sprintf $template;
    warn sprintf $template, $_ for split /\n/, $lines;
  }
  
  warn sprintf $template;
  warn "$line\n";
  warn "\n";
}

1;
