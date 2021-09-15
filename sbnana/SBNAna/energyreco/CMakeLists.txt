# Allow "StandardRecord/" to find the right directory
include_directories(..)

cet_make(LIBRARY_NAME energyreco
         LIBRARIES CAFAnaCore
                   CAFAnaSysts
                   SBNAnaVars
                   SBNAnaCuts
                   StandardRecordProxy
        )

install_headers()
install_source()

