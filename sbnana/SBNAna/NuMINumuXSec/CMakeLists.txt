# Allow "StandardRecord/" to find the right directory
include_directories( $ENV{CETLIB_INC} )
include_directories( $ENV{CETLIB_EXCEPT_INC} )
include_directories(..)

link_directories(
  $ENV{CETLIB_EXCEPT_LIB}
  $ENV{CETLIB_LIB}
)

cet_make(LIBRARY_NAME NuMINumuXSec
         LIBRARIES CAFAnaCore
                   CAFAnaSysts
                   SBNAnaVars
                   SBNAnaCuts
                   StandardRecordProxy
                   cetlib
                   cetlib_except
        )

install_headers()
install_source()
