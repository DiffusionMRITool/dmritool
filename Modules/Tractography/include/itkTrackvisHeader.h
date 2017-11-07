/**
 *       @file  itkTrackvisHeader.h
 *      @brief  
 *     Created  "07-22-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkTrackvisHeader_h
#define __itkTrackvisHeader_h

#include <iostream>
#include <string>
#include <cstring>

#include "utlCoreMacro.h"
#include "utlSmartAssert.h"

#include "utlDMRI.h"

namespace itk
{

/** 
* Trackvis header format of fiber track file. 
* http://www.trackvis.org/docs/?subsect=fileformat  
* */
typedef struct {
	char id_string[6]="TRACK";            // ID string for track file. The first 5 characters must be "TRACK".
	short int dim[3];	                    // Dimension of the image volume.
	float voxel_size[3];		              // Voxel size of the image volume.
	float origin[3];                      // Origin of the image volume. This field is not yet being used by TrackVis. That means the origin is always (0, 0, 0).
	short int n_scalars=0;                // Number of scalars saved at each track point (besides x, y and z coordinates).
  char scalar_name[10][20];             // Name of each scalar. Can not be longer than 20 characters each. Can only store up to 10 names.
  short int n_properties=0;	            // Number of properties saved at each track.
  char property_name[10][20];	          // Name of each property. Can not be longer than 20 characters each. Can only store up to 10 names.
  float vox_to_ras[4][4];	              // 4x4 matrix for voxel to RAS (crs to xyz) transformation. If vox_to_ras[3][3] is 0, it means the matrix is not recorded. This field is added from version 2. 
  char reserved[444];                   // Reserved space for future version.
  char voxel_order[4]="RAS";            // Storing order of the original image data. Explained here.
  char pad2[4];                         // Paddings.
  float image_orientation_patient[6];   // Image orientation of the original image. As defined in the DICOM header.
  char pad1[2];                         // Paddings.
  unsigned char invert_x;               // Inversion/rotation flags used to generate this track file. For internal use only.
  unsigned char invert_y;               // As above.
  unsigned char invert_z;               // As above.
  unsigned char swap_xy;                // As above.
  unsigned char swap_yz;                // As above.
  unsigned char swap_zx;                // As above.
  int n_count=0;                        // Number of tracks stored in this track file. 0 means the number was NOT stored.
  int version=2;                        // Version number. Current version is 2.
  int hdr_size=1000;                    // Size of the header. Used to determine byte swap. Should be 1000.
} TrackVisHeaderType;

void
CopyTrackvisHeader(const TrackVisHeaderType& hFrom, TrackVisHeaderType& hTo)
{
  memcpy(&hTo, &hFrom, 1000);
}

/** test if two fibers can be merged  */
bool
IsSameStructure(const TrackVisHeaderType& h1, const TrackVisHeaderType& h2)
{
  for ( int i = 0; i < 3; ++i ) 
    {
    if (h1.dim[i]!=h2.dim[i] || h1.voxel_size[i]!=h2.voxel_size[i])
      return false;
    }

  if (h1.n_scalars!=h2.n_scalars || h1.n_properties!=h2.n_properties)
    return false;

  for ( int i = 0; i < h1.n_properties; ++i ) 
    {
    if (strcmp(h1.property_name[i], h2.property_name[i])!=0)
      return false;
    }
  for ( int i = 0; i < h1.n_scalars; ++i ) 
    {
    if (strcmp(h1.scalar_name[i], h2.scalar_name[i])!=0)
      return false;
    }

  return true;
}

inline std::vector<std::string>
GetScalarNames(const TrackVisHeaderType& header)
{
  return utl::CovertChar2DArrayToStringArray(header.scalar_name, 10);
}

inline std::vector<std::string>
GetPropertyNames(const TrackVisHeaderType& header)
{
  return utl::CovertChar2DArrayToStringArray(header.property_name, 10);
}

inline bool
HasPropertyName(const TrackVisHeaderType& header, const std::string& name)
{
  for ( int i = 0; i < header.n_properties; ++i ) 
    {
    if (name==std::string(header.property_name[i]))
      return true;
    }
  return false;
}

inline bool
HasScalarName(const TrackVisHeaderType& header, const std::string& name)
{
  for ( int i = 0; i < header.n_scalars; ++i ) 
    {
    if (name==std::string(header.scalar_name[i]))
      return true;
    }
  return false;
}

/** Get total dimension of scalars  */
inline int
GetDimensionOfScalars(const TrackVisHeaderType& header)
{
  int num=0;
  for ( int i = 0; i < header.n_scalars; ++i ) 
    {
    std::string name(header.scalar_name[i]);
    if (name=="")
      return num;
    else
      num += utl::GetScalarsDimentionByName(name);
    }
  return num;
}

/** Get total dimension of properties  */
inline int
GetDimensionOfProperties(const TrackVisHeaderType& header)
{
  int num=0;
  for ( int i = 0; i < header.n_properties; ++i ) 
    {
    std::string name(header.property_name[i]);
    if (name=="")
      return num;
    else
      num += utl::GetScalarsDimentionByName(name);
    }
  return num;
}

/** remove scalars by its name   */
inline void
RemoveScalarName(TrackVisHeaderType& header, const std::string& name)
{
  int i=0;
  for ( i = 0; i < header.n_scalars; ++i ) 
    {
    if (name==std::string(header.scalar_name[i]))
      break;
    }

  for ( int j = i; j+1 < header.n_scalars; ++j ) 
    strcpy(header.scalar_name[j], header.scalar_name[j+1]);
  strcpy(header.scalar_name[header.n_scalars-1], "");

  header.n_scalars--;
}

/** remove properties by its name   */
inline void
RemovePropertyName(TrackVisHeaderType& header, const std::string& name)
{
  int i=0;
  for ( i = 0; i < header.n_properties; ++i ) 
    {
    if (name==std::string(header.property_name[i]))
      break;
    }

  for ( int j = i; j+1 < header.n_properties; ++j ) 
    strcpy(header.property_name[j], header.property_name[j+1]);
  strcpy(header.property_name[header.n_properties-1], "");

  header.n_properties--;
}


/** read trackvis header  */
inline void
ReadTrackVisHeader( const std::string& filename, TrackVisHeaderType& header )
{
  // read header
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  if (!file)
    utlGlobalException(true, "Unable to open file " + filename);
  
  // read header
  fseek (file, 0, SEEK_SET);  
  int hsize =  fread((char*)(&header), 1, 1000, file);
  utlGlobalException(hsize!=1000, "reading ERROR");
  utlGlobalException(strcmp(header.id_string,"TRACK")!=0, "Wrong data format. Magic code is wrong. It is " + std::string(header.id_string) +". Should be 'TRACK'");

  // In case n_count is not correct, read remaining data file to get the correct n_count
  fseek (file, 0, SEEK_END);  
  long fsize = ftell(file);
  long offset=1000;
  int n_count=0; 
  int dim_scalars=GetDimensionOfScalars(header);
  int dim_properties=GetDimensionOfProperties(header);
  while (offset<fsize)
    {
    fseek (file, offset, SEEK_SET);
    int numPoints;
    fread((char*)&numPoints, sizeof(int), 1, file);
    offset += 4+ numPoints*(3+dim_scalars)*4 + dim_properties*4;
    n_count++;
    }

  // If n_count is not correct in header, then set it as the correct value from the file.
  if (n_count!=header.n_count)
    {
    utlSAGlobalWarning(n_count!=header.n_count)(n_count)(header.n_count).
      msg("The number of tracts in the header does not match the actual number of tracts in the file. We set n_count as the actual number of tracts in the file");
    header.n_count = n_count;
    }
  fclose (file);
}

/** write trackvis header  */
inline void
WriteTrackVisHeader( const TrackVisHeaderType& header, FILE* file )
{
  fseek(file, 0, SEEK_SET);
  int hsize = fwrite((char*)&header, 1, 1000, file);
  utlGlobalException (hsize!=1000, "Error when saving the fiber!");
}



/** print trackvis header  */
inline void 
PrintTractVisHeader( const TrackVisHeaderType& header, std::ostream & os=std::cout) 
{
  os << "id_string                 : " << header.id_string  << std::endl;
  os << "version                   : " << header.version  << std::endl;
  os << "dim                       : " << header.dim[0] << ", " << header.dim[1] << ", " << header.dim[2] << std::endl;
  os << "voxel_size                : " << header.voxel_size[0] << ", " << header.voxel_size[1] << ", " << header.voxel_size[2] << std::endl;
  os << "origin                    : " << header.origin[0] << ", " << header.origin[1] << ", " << header.origin[2] << std::endl;
  os << "number of tracks          : " << header.n_count << std::endl;
  os << "number of scalars         : " << header.n_scalars  << std::endl;
  for ( int i = 0; i < header.n_scalars; ++i ) 
    {
      os << "scalar_name[" << i << "]            : " << std::string(header.scalar_name[i]) << std::endl;
    }
  os << "number of properties      : " << header.n_properties  << std::endl;
  os << "dimension of scalars      : " << GetDimensionOfScalars(header)  << std::endl;
  os << "dimension of properties   : " << GetDimensionOfProperties(header)  << std::endl;
  for ( int i = 0; i < header.n_properties; ++i ) 
    {
      os << "property_name[" << i << "]          : " << std::string(header.property_name[i]) << std::endl;
    }
  if (header.vox_to_ras[3][3]!=0)
    {
    os << "vox_to_ras                : " << header.vox_to_ras[0][0] << ", " << header.vox_to_ras[0][1] << ", " << header.vox_to_ras[0][2] << ", " << header.vox_to_ras[0][3] << std::endl;
    os << "                            " << header.vox_to_ras[1][0] << ", " << header.vox_to_ras[1][1] << ", " << header.vox_to_ras[1][2] << ", " << header.vox_to_ras[1][3] << std::endl;
    os << "                            " << header.vox_to_ras[2][0] << ", " << header.vox_to_ras[2][1] << ", " << header.vox_to_ras[2][2] << ", " << header.vox_to_ras[2][3] << std::endl;
    os << "                            " << header.vox_to_ras[3][0] << ", " << header.vox_to_ras[3][1] << ", " << header.vox_to_ras[3][2] << ", " << header.vox_to_ras[3][3] << std::endl;
    }
  os << "voxel_order               : " << header.voxel_order  << std::endl;

  const float* a= header.image_orientation_patient;
  os << "image_orientation         : " << a[0]<<","<<a[1]<<","<<a[2]<<","<<a[3]<<","<<a[4]<<","<<a[5]<<   std::endl;
}

}


#endif 
