#pragma once

#include <memory>
#include <map>
#include <cmath>

#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/vec2.hpp>

#include "ao/ui/gl/core.hpp"
#include "ao/ui/gl/frame.hpp"
#include "ao/ui/gl/axes.hpp"
#include "ao/ui/gl/border.hpp"

class Tree;

class Window
{
public:
    /*
     *  We should only have one window constructed at a time
     */
    static Window* singleton;
    static Window* instance();

    /*
     *  Window destructor deletes all associated Frames
     */
    ~Window();

    /*
     *  Adds a Frame to redraw the given shape asynchronously
     *  (will not take effect until the next poll() call)
     */
    void addTree(std::string filename, std::string name, Tree* t);

    /*
     *  Clears existing frames asynchronously
     *  (will not take effect until the next poll() call)
     */
    void clearFrames();

    /*
     *  Clears frames associated with the given filename
     *  Executed synchronously (so must be called from main thread)
     */
    void clearFile(std::string filename);

    /*
     *  Runs the window forever
     */
    void run();

    /*
     *  Polls for window events
     */
    void poll();

    /*
     *  Redraw the window with the current state
     */
    void draw() const;

protected:
    /*
     *  Default Window constructor
     *
     *  Triggers an assertion failure on OpenGL errors
     */
    explicit Window();

    /*
     *  Check if the incoming pointer is populated and create a Frame if so
     *
     *  Returns true if the window should be redrawn.
     */
    bool loadFrame();

    /*
     *  Iterate over all frames, checking their render status.
     *
     *  Returns true if the window should be redrawn.
     */
    bool pollFrames();

    /*
     *  Checks if the clear flag is set, clearing frames if it is
     *  (by moving them to the stale list)
     *
     *  Returns true if the window should be redrawn.
     */
    bool checkClear();

    /*
     *  Checks the frames in the stale list, deleting those that are done
     */
    void pruneStale();

    /*
     *  Request that every child Frame render itself
     */
    void render();

    /*
     *  Projection matrix
     *  (compensates for window size)
     */
    glm::mat4 proj() const;

    /*
     *  View matrix
     *  Applies scale, rotation, offset
     */
    glm::mat4 view() const;

    /*
     *  Complete transform matrix
     */
    glm::mat4 M() const { return proj() * view(); }

    /*
     *  Returns a mouse position scaled to the range -1, 1
     */
    glm::vec2 scaledMousePos(double x, double y) const;

    /*
     *  Returns a mouse position in integer screen coordinates
     */
    glm::ivec2 rawMousePos(double x, double y) const;

    /*
     *  Resized window callback
     */
    static void _resized(GLFWwindow* window, int w, int h);
    void resized(int w, int h);

    /*
     *  Callback for mouse movement
     */
    static void _mouseMove(GLFWwindow* window, double x, double y);
    void mouseMove(double x, double y);

    /*
     *  Callback for mouse scroll
     */
    static void _mouseButton(GLFWwindow* window, int b, int a, int m);
    void mouseScroll(double x, double y);

    /*
     *  Mouse press callback
     */
    static void _mouseScroll(GLFWwindow* window, double sx, double sy);
    void mouseButton(int button, int action, int mods);

    /*  Pointer to raw window  */
    GLFWwindow* const window;

    /*  Atomic pointer to incoming tree  */
    std::atomic<std::tuple<std::string, std::string, Tree*>*> incoming;

    /*  Set to true if frames should be cleared  */
    std::atomic_bool clear;

    /*  Width and height are of the framebuffer, rather than the window  *
     *  (to properly cope with high DPI monitors)                        */
    int width;
    int height;

    /*  Store the scene's center, scale, and rotation                    *
     *  (Euler angles aren't the best representation, but they're easy)  */
    glm::vec3 center;
    float scale=1;
    float pitch=0;
    float  roll=0;

    enum { DRAG_NONE,
           DRAG_PAN,
           DRAG_ROTATE,
           DRAG_MOVE_WINDOW,
           DRAG_SCALE_WINDOW} drag_mode=DRAG_NONE;

    /*  This is the current mouse position (in -1, 1 window coordinates)  */
    glm::vec2 mouse_pos;

    /*  This is the current mouse position in global coordinates    *
     *  (where 0,0 is the top left corner of the screen)            */
    glm::ivec2 cursor_pos;

    /*  Objects to draw in 3D viewport  */
    Axes axes;
    Border border;
    std::map<std::string, std::map<std::string, Frame*>> frames;

    /*  Shader strings to draw window border  */
    std::string vert;
    std::string frag;

    /*  Frames that should be deleted once they are done rendering  */
    std::list<Frame*> stale;

    /*  Border for dragging window around  */
    static const int BORDER=30;
};
