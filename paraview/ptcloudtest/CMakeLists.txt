ttk_add_paraview_plugin(ttkptcloudtest
	SOURCES ${VTKWRAPPER_DIR}/ttkptcloudtest/ttkptcloudtest.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML ptcloudtest.xml
    LINK common pointDistField)

