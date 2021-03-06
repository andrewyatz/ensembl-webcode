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

package EnsEMBL::Web::ExtIndex::PFETCH;

use strict;
use IO::Socket;
use Sys::Hostname;

sub new { my $class = shift; my $self = {}; bless $self, $class; return $self; }

sub get_seq_by_acc { my $self = shift; $self->get_seq_by_id( @_ ); }

sub get_seq_by_id {
  my ($self, $arghashref)=@_;

  # Get the ID to pfetch
  my $str = $arghashref->{ID} || return [];

  # Additional options
  if( $arghashref->{OPTIONS} eq 'desc'       ) { $str .= " -D" }
  if( $arghashref->{OPTIONS} =~ /(-d\s+\w+)/ ) { $str .= " $1" }
  if( $arghashref->{DB} eq 'PUBLIC'          ) { $str .= " -d public" }
  if( $arghashref->{DB} =~ /UNIPROT/         ) { $str = " -a $str" }

  # Get the pfetch server
  my $server = $self->fetch_pfetch_server(
    $arghashref->{'species_defs'}->ENSEMBL_PFETCH_SERVER,
    $arghashref->{'species_defs'}->ENSEMBL_PFETCH_PORT
  );

  my $hostname = &Sys::Hostname::hostname();
  my $output;
  if ($arghashref->{'strand_mismatch'}) {
    print $server "--client $hostname $str -r \n";
    push @$output, $_ while(<$server>);
  }
  else {
    print $server "--client $hostname $str  \n";
    push @$output, $_ while(<$server>);
  }

  return $output;
}

sub fetch_pfetch_server {
  my $self   = shift;
  my $server = shift;
  my $port   = shift;
  if( ! $server ){ die "No ENSEMBL_PFETCH_SERVER found in config" }
  
  my $s = IO::Socket::INET->new( PeerAddr => $server,
    PeerPort => $port, Proto    => 'tcp', Type     => SOCK_STREAM, Timeout  => 10,
  );
  if ($s){
    $s->autoflush(1);
    return( $s );
  } 
  die "Cannot connect to the Trace server - please try again later.";
}

1;
