
import vtk


def convertDataSetToSurface(algorithmOutputPort):
    dataSetSurfaceFilter = vtk.vtkDataSetSurfaceFilter()
    dataSetSurfaceFilter.SetInputConnection(algorithmOutputPort)
    dataSetSurfaceFilter.UseStripsOn()
    dataSetSurfaceFilter.Update()
    polyData = vtk.vtkPolyData()
    polyData.ShallowCopy(dataSetSurfaceFilter.GetOutput())
    return polyData


def readLegacyVTK(file_name):
    '''Support POLYDATA and some other vtk data (*.vtk).'''
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(file_name)
    reader.Update()
    if None != reader.GetPolyDataOutput():
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(reader.GetPolyDataOutput())
        return polyData, 'PolyData'
    if None != reader.GetUnstructuredGridOutput():
        return convertDataSetToSurface(reader.GetOutputPort()), 'UnstructuredGrid'
    if None != reader.GetStructuredPointsOutput():
        return convertDataSetToSurface(reader.GetOutputPort()), 'StructuredPoints'
    if None != reader.GetStructuredGridOutput():
        return convertDataSetToSurface(reader.GetOutputPort()), 'StructuredGrid'
    if None != reader.GetRectilinearGridOutput():
        return convertDataSetToSurface(reader.GetOutputPort()), 'RectilinearGrid'
    else:
        raise Exception("Unsupported type!\n")


def readDataSet(file_name, readerType):
    reader = readerType()
    reader.SetFileName(file_name)
    reader.Update()
    return convertDataSetToSurface(reader.GetOutputPort())


def readPolydata(file_name):
    """ Load a vtk polydata to a supported format file

    Parameters
    ----------
        file_name : string

    Returns
    -------
        output : vtkPolyData
    """
    # get file extension (type) lower case
    file_extension = file_name.split(".")[-1].lower()

    if file_extension == "vtk":
        try:
            return readLegacyVTK(file_name)[0]
        except:
            reader = vtk.vtkPolyDataReader()
    elif file_extension == "fib":
        reader = vtk.vtkPolyDataReader()
    elif file_extension == "ply":
        reader = vtk.vtkPLYReader()
    elif file_extension == "stl":
        reader = vtk.vtkSTLReader()
    elif file_extension == "vtp":
        reader = vtk.vtkXMLPolyDataReader()
    elif file_extension == "vti":
        return readDataSet(file_name, vtk.vtkXMLImageDataReader)
    elif file_extension == "vts":
        return readDataSet(file_name, vtk.vtkXMLStructuredGridReader)
    elif file_extension == "vtu":
        return readDataSet(file_name, vtk.vtkXMLUnstructuredGridReader)
    elif file_extension == "vtr":
        return readDataSet(file_name, vtk.vtkXMLRectilinearGridReader)
    elif file_extension == "xml":
        reader = vtk.vtkXMLPolyDataReader()
    elif file_extension == "obj":
        try:  # try to read as a normal obj
            reader = vtk.vtkOBJReader()
        except:  # than try load a MNI obj format
            reader = vtk.vtkMNIObjectReader()
    else:
        raise Exception("polydata " + file_extension + " is not suported")

    reader.SetFileName(file_name)
    reader.Update()
    return reader.GetOutput()


def savePolydata(polydata, file_name, binary=False, color_array_name=None):
    """ Save a vtk polydata to a supported format file

    Parameters
    ----------
        polydata : vtkPolyData
        file_name : string
    """
    # get file extension (type)
    file_extension = file_name.split(".")[-1].lower()

    if file_extension == "vtk":
        writer = vtk.vtkPolyDataWriter()
    elif file_extension == "fib":
        writer = vtk.vtkPolyDataWriter()
    elif file_extension == "ply":
        writer = vtk.vtkPLYWriter()
    elif file_extension == "stl":
        writer = vtk.vtkSTLWriter()
    elif file_extension == "xml":
        writer = vtk.vtkXMLPolyDataWriter()
    elif file_extension == "vtp":
        writer = vtk.vtkXMLPolyDataWriter()
    elif file_extension == "vti":
        writer = vtk.vtkXMLImageDataWriter()
    elif file_extension == "vts":
        writer = vtk.vtkXMLStructuredGridWriter()
    elif file_extension == "vtu":
        writer = vtk.vtkXMLUnstructuredGridWriter()
    elif file_extension == "vtr":
        writer = vtk.vtkXMLRectilinearGridWriter()
    elif file_extension == "obj":
        raise Exception("mni obj or Wavefront obj ?")

    writer.SetFileName(file_name)
    writer.SetInputData(polydata)
    if color_array_name is not None:
        writer.SetArrayName(color_array_name)

    if binary:
        writer.SetFileTypeToBinary()
    writer.Update()
    writer.Write()


def flipVTKImageData(image, image_flip):
    '''flip image for visualization.

    Parameters
    ----------
        image : vtkImageData
        image_flip : tuple, (-1 1 1) means flip x-axis

    Returns
    -------
        imageOut : vtkImageData, flipped image
    '''

    if image_flip[0] == 1 and image_flip[1] == 1 and image_flip[2] == 1:
        return image

    imageOut = vtk.vtkImageData()
    imageOut.DeepCopy(image)

    for i, ff in enumerate(image_flip):
        if ff == -1:
            flipFilter = vtk.vtkImageFlip()
            flipFilter.SetInputData(imageOut)
            flipFilter.SetFilteredAxis(i)
            flipFilter.Update()
            imageOut = flipFilter.GetOutput()

    return imageOut


def printImage(image):
    '''print image'''
    dim = image.GetDimensions()
    for i in range(dim[0]):
        for j in range(dim[1]):
            for k in range(dim[2]):
                if image.GetScalarComponentAsFloat(i, j, k, 0) > 0:
                    print i, j, k, image.GetScalarComponentAsFloat(i, j, k, 0)

