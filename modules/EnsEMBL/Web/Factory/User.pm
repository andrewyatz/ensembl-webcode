package EnsEMBL::Web::Factory::User;

use strict;

use EnsEMBL::Web::Factory;
use EnsEMBL::Web::DBSQL::UserDB;
use EnsEMBL::Web::ParameterSet;
use EnsEMBL::Web::Proxy::Object;
our @ISA = qw(EnsEMBL::Web::Factory);

sub createObjects { 
  my $self          = shift;
  warn "FACTORY new USER: " . $ENV{'ENSEMBL_USER_ID'};
  my $adaptor = EnsEMBL::Web::DBSQL::UserDB->new();
  my $parameter_set = EnsEMBL::Web::ParameterSet->new((
                                      cgi => $self->__data->{'_input'},
                                                )); 

  $self->DataObjects( new EnsEMBL::Web::Proxy::Object
    (
     'User', { adaptor => $adaptor, id => $ENV{'ENSEMBL_USER_ID'}, parameter_set => $parameter_set }, $self->__data
    ) 
  );

}

1;
