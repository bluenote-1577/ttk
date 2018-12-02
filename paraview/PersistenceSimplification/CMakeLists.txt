ttk_add_paraview_plugin(ttkPersistenceSimplification
	SOURCES ${VTKWRAPPER_DIR}/ttkPersistenceSimplification/ttkPersistenceSimplification.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML PersistenceSimplification.xml
	LINK persistenceSimplification)

