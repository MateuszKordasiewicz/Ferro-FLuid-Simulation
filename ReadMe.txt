/////////////////////////////////////////////////////////////////////////////////////
//             Project Code 'code/'
/////////////////////////////////////////////////////////////////////////////////////
Project code is contained in the code/ folder.
This contains the geometry folder which has the .geo file for the mesh and the shell script for running the mesh maker.
This also contains an output folder which is used when the code is ran without the mesh or time study.
Then the mesh and time study folders are used if the respective booleans are set to 'true'.

/////////////////////////////////////////////////////////////////////////////////////
//             Paraview animations 'movies/'
/////////////////////////////////////////////////////////////////////////////////////
This folder contains two sub folders, one with movies from the mesh study and another with movies for the results section of the report.
The results section folder also has two sub folders, one for the discussion of the finest mesh and another for the second discussion about larger timesteps.


/////////////////////////////////////////////////////////////////////////////////////
//             Output analysis 'python_code/'
/////////////////////////////////////////////////////////////////////////////////////
This contains two folders one for the time study and one for the mesh study.
In both of these folders the data outputs from the mesh and time studies is contained along with the python code used to do the analysis.
There is also a VTKs folder which contains all of the error vectors exported as VTKs for viewing.