/*
Studio: a simple GUI for the libfive CAD kernel
Copyright (C) 2017  Matt Keeter

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/
#include <QColor>
#include <Eigen/StdVector>
#include <boost/math/constants/constants.hpp>

#include "studio/arrow.hpp"
#include "studio/shader.hpp"

namespace Studio {

Arrow::Arrow()
{
    // Nothing to do here
}

void Arrow::initializeGL(int res)
{
    initializeOpenGLFunctions();

    tri_count = res * 8;
    Eigen::Array<float, Eigen::Dynamic, 18, Eigen::RowMajor> vs(tri_count, 18);

    Eigen::Array<float, 1, 18> zflip = Eigen::Array<float, 1, 18>::Ones();
    for (unsigned i=2; i < 18; i += 3)
    {
        zflip(i - 1) = -1;
        zflip(i) = -1;
    }

    const float pi = boost::math::constants::pi<float>();
    const float ro = 0.2;
    const float ri = 0.1;
    const float ah = 0.4;
    for (int i=0; i < res; ++i)
    {
        float a0 = pi * 2 * i / (res - 1.0);
        float a1 = pi * 2 * (i + 1) / (res - 1.0);
        float x0 = cos(a0);
        float y0 = sin(a0);
        float x1 = cos(a1);
        float y1 = sin(a1);

        // Each row is a triangle, with per-vertex values for
        //                 Position                 Normal
        vs.row(8*i + 0) << 0,0,1,                   (x0 + x1),(y0 + y1),0,
                           x0*ro, y0*ro, 1 - ah,    x0, y0, 0,
                           x1*ro, y1*ro, 1 - ah,    x1, y1, 0;
        vs.row(8*i + 1) << x0*ro, y0*ro, 1 - ah,    0, 0, -1,
                           x1*ri, y1*ri, 1 - ah,    0, 0, -1,
                           x1*ro, y1*ro, 1 - ah,    0, 0, -1;
        vs.row(8*i + 2) << x0*ro, y0*ro, 1 - ah,    0, 0, -1,
                           x0*ri, y0*ri, 1 - ah,    0, 0, -1,
                           x1*ri, y1*ri, 1 - ah,    0, 0, -1;
        vs.row(8*i + 3) << x1*ri, y1*ri, 1 - ah,    x1, y1, 0,
                           x0*ri, y0*ri, 1 - ah,    x0, y0, 0,
                           x0*ri, y0*ri, ah - 1,    x0, y0, 0;
        vs.row(8*i + 4) << x1*ri, y1*ri, 1 - ah,    x1, y1, 0,
                           x0*ri, y0*ri, ah - 1,    x0, y0, 0,
                           x1*ri, y1*ri, ah - 1,    x0, y0, 0;

        vs.row(8*i + 5) = vs.row(8*i + 0) * zflip;
        vs.row(8*i + 6) = vs.row(8*i + 1) * zflip;
        vs.row(8*i + 7) = vs.row(8*i + 2) * zflip;
    }

    vao.create();
    vao.bind();

    vbo.create();
    vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
    vbo.bind();
    vbo.allocate(vs.data(), vs.size() * sizeof(*vs.data()));

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(GLfloat), NULL);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 6*sizeof(GLfloat),
                          (GLvoid*)(3 * sizeof(GLfloat)));
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(2);

    vbo.release();
    vao.release();
}

void Arrow::draw(QMatrix4x4 M, QVector3D pos, float scale,
                 QVector3D norm, QColor color)
{
    // Build an orthonormal basis, per Duff 2017
    QMatrix4x4 rot;
    rot.setColumn(2, norm.normalized().toVector4D());
    if (norm.z() < 0.)
    {
        const float a = 1.0f / (1.0f - norm.z());
        const float b = norm.x() * norm.y() * a;
        rot.setColumn(0, {1.0f - norm.x() * norm.x() * a, -b, norm.x(), 0});
        rot.setColumn(1, {b, norm.y() * norm.y()*a - 1.0f, -norm.y(), 0});
    }
    else
    {
        const float a = 1.0f / (1.0f + norm.z());
        const float b = -norm.x() * norm.y() * a;
        rot.setColumn(0, {1.0f - norm.x() * norm.x() * a, b, -norm.x(), 0});
        rot.setColumn(1, {b, 1.0f - norm.y() * norm.y() * a, -norm.y(), 0});
    }

    M.translate(pos);
    M.scale(scale);
    M *= rot;

    Shader::basic->bind();
    glUniformMatrix4fv(Shader::basic->uniformLocation("M"),
                       1, GL_FALSE, M.data());
    glUniform1i(Shader::basic->uniformLocation("shading"), 1); // per-triangle
    glUniform4f(Shader::basic->uniformLocation("color_add"),
                color.redF()*0.6, color.greenF()*0.6, color.blueF()*0.6, 0.0f);
    vao.bind();

    // Draw once at 20% opacity but without depth culling
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_STENCIL_TEST);
    glStencilFuncSeparate(GL_FRONT, GL_GREATER, 1, 0xFF);
    glStencilFuncSeparate(GL_BACK, GL_NEVER, 0, 0xFF);
    glStencilOpSeparate(GL_FRONT, GL_KEEP, GL_INCR, GL_INCR);
    glStencilOpSeparate(GL_BACK, GL_KEEP, GL_INCR, GL_INCR);

    glUniform4f(Shader::basic->uniformLocation("color_mul"),
                color.redF()/2, color.greenF()/2, color.blueF()/2, 0.3f);
    glDrawArrays(GL_TRIANGLES, 0, tri_count * 3);

    // Then draw in full color, with depth culling on
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_STENCIL_TEST);

    glUniform4f(Shader::basic->uniformLocation("color_mul"),
                color.redF()/2, color.greenF()/2, color.blueF()/2, 1.0f);
    glDrawArrays(GL_TRIANGLES, 0, tri_count * 3);

    vao.release();
    Shader::basic->release();
}

}   // namespace Studio
