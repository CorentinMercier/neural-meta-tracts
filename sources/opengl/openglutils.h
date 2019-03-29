#ifndef OPENGLUTILS_H
#define OPENGLUTILS_H

#include <GL/glew.h>
#include <GL/gl.h>
#include <string>
#include <iostream>
#include <fstream>
#include <streambuf>

using namespace std;

/** Macro indicating OpenGL errors (file and line) */
#define PRINT_OPENGL_ERROR() print_opengl_error(__FILE__, __LINE__)

typedef struct
{
    GLenum type;
    const char* filename;
    GLuint shader;
}ShaderInfo;

//Function returning the ID of the program using the shaders included in ShaderInfo
GLuint loadShaders(ShaderInfo* shaders);

/** Print OpenGL Error given the information of the file and the line
 *  Function called by the macro PRINT_OPENGL_ERROR */
bool print_opengl_error(const char *file, int line);

GLint get_uni_loc(GLuint program, const GLchar *name);


#endif // OPENGLUTILS_H
