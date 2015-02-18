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

#include <QObject>
#include <QString>
#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QShortcut>

#include "VTKViewer.h"

static void makeShortcut(
    int key, QWidget & parent, QWidget & target, const char * slot) {
  QObject::connect(
      new QShortcut(QKeySequence(key), &parent),
      SIGNAL(activated()), &target, slot);
}

int main(int argc, char ** argv) {
  if (argc == 1)
    {
    std::cerr
      << "\nUseage:  \n  " << argv[0]
      << " FILE [MORE FILES...]\n"
      "Supported File Formats:\n"
      "  *.vtk - VTK Legacy File\n"
      "  *.vtp - VTK Polygonal Data File\n"
      "  *.ply - Stanford Polygon File\n"
      "  *.obj - Wavefront Object file\n"
      "  *.stl - Stereolithography File\n"
      "  *.pdb - Protein Data Bank File\n"
      "  *.vtu - VTK Unstructured Grid Data File\n"
      "Controls:\n"
      "  's' - surface\n"
      "  'w' - wireframe\n"
      "  'r' - reset and center camera\n"
      "  'ctrl-q' - quit\n"
      "  'ctrl-r' - toggle rotation\n"
      "  'ctrl-s' - toggle stereo mode\n"
      "  'ctrl-t' - change stereo type\n"
      "  'ctrl-p' - screenshot\n"
      "More Info:\n"
      "  https://github.com/HalCanary/vtkviewer\n\n";
    return 1;
    }

  QApplication app(argc, argv);
  QMainWindow mw;
  mw.setWindowTitle("VTK Viewer");
  VTKViewer v;

  makeShortcut(Qt::CTRL + Qt::Key_Q, mw, mw, SLOT(close()));
  makeShortcut(Qt::CTRL + Qt::Key_R, mw, v, SLOT(toggleRotate()));
  makeShortcut(Qt::CTRL + Qt::Key_S, mw, v, SLOT(toggleStereo()));
  makeShortcut(Qt::CTRL + Qt::Key_T, mw, v, SLOT(nextStereoType()));
  makeShortcut(Qt::CTRL + Qt::Key_P, mw, v, SLOT(screenshot()));

  mw.setCentralWidget(&v);
  for (int i = 1; i < argc; ++i)
    {
    v.add(argv[i]);
    }
  // mw.showMaximized();
  mw.show();
  return app.exec();
}
