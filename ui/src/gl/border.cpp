#include <cassert>

#include "ao/ui/gl/border.hpp"
#include "ao/ui/gl/shader.hpp"

// Vertex shader
const std::string Border::vert = R"(
#version 330

layout(location=0) in vec3 vertex_position;

void main()
{
    gl_Position = vec4(vertex_position, 1.0f);
}
)";

// Fragment shader
const std::string Border::frag = R"(
#version 330

uniform int w;
uniform int h;
uniform int b;

out vec4 fragColor;

void main()
{
    vec4 c = vec4(0.2f, 0.2f, 0.2f, 1.0f);

    if (gl_FragCoord.y > h - b)
    {
        fragColor = c;
    }
    else if (gl_FragCoord.x > w - b && gl_FragCoord.y < b &&
             gl_FragCoord.x - w + b > gl_FragCoord.y)
    {
        fragColor = c;
    }
    else
    {
        discard;
    }
}
)";

////////////////////////////////////////////////////////////////////////////////

Border::Border()
    : vs(Shader::compile(vert, GL_VERTEX_SHADER)),
      fs(Shader::compile(frag, GL_FRAGMENT_SHADER)),
      prog(Shader::link(vs, fs))
{
    assert(vs);
    assert(fs);
    assert(prog);

    glGenBuffers(1, &vbo);
    glGenVertexArrays(1, &vao);

    glBindVertexArray(vao);
    {
        GLfloat vertices[] = {-1.0f, -1.0f, 0.0f,
                               1.0f, -1.0f, 0.0f,
                               1.0f,  1.0f, 0.0f,
                              -1.0f,  1.0f, 0.0f};
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices),
                     vertices, GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE,
                              3 * sizeof(GLfloat), (GLvoid*)0);
        glEnableVertexAttribArray(0);
    }
    glBindVertexArray(0);
}

void Border::draw(int w, int h, int border) const
{
    // Active the shader
    glUseProgram(prog);

    // Bind the vertex array (which sets up VBO bindings)
    glBindVertexArray(vao);

    glUniform1i(glGetUniformLocation(prog, "w"), w);
    glUniform1i(glGetUniformLocation(prog, "h"), h);
    glUniform1i(glGetUniformLocation(prog, "b"), border);

    // Draw the quad!
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);

    glBindVertexArray(0);
}
