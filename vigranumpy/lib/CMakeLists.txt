SET(SOURCES
    __init__.py
    arraytypes.py
    ufunc.py)

INSTALL(FILES ${SOURCES} DESTINATION ${VIGRANUMPY_INSTALL_DIR}/vigra)
