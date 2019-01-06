ttk_add_vtk_library(ttkScalarFieldFromPointCloud
    SOURCES ttkScalarFieldFromPointCloud.cpp
	HEADERS ttkScalarFieldFromPointCloud.h
    LINK common scalarFieldFromPointCloud ttkTriangulation)
