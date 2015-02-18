// VTK Viewer
// Written 2012 Hal Canary <http://cs.unc.edu/~hal>
// Copyright 2012 University of North Carolina at Chapel Hill.
//
// Licensed under the Apache License, Version 2.0 (the "License"); you
// may not use this file except in compliance with the License.  You
// may obtain a copy of the License at
//
//   LICENSE.md in this repository or
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
// implied.  See the License for the specific language governing
// permissions and limitations under the License.
#ifndef VTKVIEWER_H
#define VTKVIEWER_H
#include <cassert>
#include <QTimer>
#include <QVTKWidget.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
class vtkPolyData;

class VTKViewer : public QVTKWidget {
  Q_OBJECT;
public:
  VTKViewer();
  void add(vtkPolyData * polyData);
  void add(const char * file_name);
public slots:
  void rotate();
  void toggleRotate();
  void toggleStereo();
  void setStereoType(int vtkStereoType);
  void nextStereoType();
  void screenshot();
  void setsize(double s1, double s2);
private:
  QTimer m_timer;
  vtkSmartPointer < vtkRenderer > m_renderer;

  double m_size[2];
};

#endif /* VTKVIEWER_H */
