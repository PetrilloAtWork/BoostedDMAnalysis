
art_make( 
          MODULE_LIBRARIES larcorealg_Geometry
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
                        canvas
			${MF_MESSAGELOGGER}
	                ${FHICLCPP}
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
        )

# install_headers()
install_fhicl()
install_source()

