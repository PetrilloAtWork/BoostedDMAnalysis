# this will create both the library and the dictionary
art_make(
  LIB_LIBRARIES
    lardataobj_RecoBase
    ${ROOT_EG} # TDatabasePDG
    ${ROOT_GENVECTORS}
    ${ROOT_CORE}
  DICT_LIBRARIES
    nusimdata_SimulationBase
  )

install_source()
install_headers()
