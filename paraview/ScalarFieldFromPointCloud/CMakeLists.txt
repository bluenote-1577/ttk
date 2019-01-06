ttk_add_paraview_plugin(ttkScalarFieldFromPointCloud
    SOURCES ${VTKWRAPPER_DIR}/ttkScalarFieldFromPointCloud/ttkScalarFieldFromPointCloud.cpp 
	PLUGIN_VERSION 1.0
    PLUGIN_XML ScalarFieldFromPointCloud.xml
    LINK common scalarFieldFromPointCloud)

