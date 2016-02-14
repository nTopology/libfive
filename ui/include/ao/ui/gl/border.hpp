#pragma once

#include "ao/ui/gl/core.hpp"

/*
 *  The Border class draws a window border
 */
class Border
{
public:
    Border();
    void draw(int w, int h, int border) const;

protected:
    // Shader strings
    static const std::string vert;
    static const std::string frag;

    GLuint vs;  // Vertex shader
    GLuint fs;  // Fragment shader
    GLuint prog;    // Shader program

    GLuint vbo; // Vertex buffer object
    GLuint vao; // Vertex array object
};
