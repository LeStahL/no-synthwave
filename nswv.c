/* No Synthwave - PC4k Intro by Team210 at Deadline 2k19
 * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "Windows.h"
#include "GL/GL.h"
#include "glext.h"

PFNGLCREATESHADERPROC glCreateShader;
PFNGLCREATEPROGRAMPROC glCreateProgram;
PFNGLSHADERSOURCEPROC glShaderSource;
PFNGLCOMPILESHADERPROC glCompileShader;
PFNGLATTACHSHADERPROC glAttachShader;
PFNGLLINKPROGRAMPROC glLinkProgram;
PFNGLUSEPROGRAMPROC glUseProgram;
PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation;
PFNGLUNIFORM1FPROC glUniform1f;

size_t strlen(const char *str)
{
	int len = 0;
	while(str[len] != '\0') ++len;
	return len;
}

void *memset(void *ptr, int value, size_t num)
{
	for(int i=num-1; i>=0; i--)
		((unsigned char *)ptr)[i] = value;
	return ptr;
}

int WINAPI demo(HINSTANCE hInstance, HINSTANCE hPrevInstance, PWSTR pCmdLine, int nCmdShow)
{
    // Display demo window
	CHAR WindowClass[]  = "Team210 Demo Window";

	WNDCLASSEX wc = { 0 };
	wc.cbSize = sizeof(wc);
	wc.style = CS_OWNDC | CS_VREDRAW | CS_HREDRAW;
	wc.lpfnWndProc = &DefWindowProc;
	wc.cbClsExtra = 0;
	wc.cbWndExtra = 0;
	wc.hInstance = hInstance;
	wc.hIcon = LoadIcon(NULL, IDI_WINLOGO);
	wc.hCursor = LoadCursor(NULL, IDC_ARROW);
	wc.hbrBackground = NULL;
	wc.lpszMenuName = NULL;
	wc.lpszClassName = WindowClass;
	wc.hIconSm = NULL;

	RegisterClassEx(&wc);
    
    HWND hwnd = CreateWindowEx(0, WindowClass, ":: Team210 :: GO - MAKE A DEMO ::", WS_POPUP | WS_VISIBLE, 0, 0, 1280, 720, NULL, NULL, hInstance, 0);
    
    DEVMODE dm = { 0 };
    dm.dmSize = sizeof(dm);
    dm.dmPelsWidth = 1280;
    dm.dmPelsHeight = 720;
    dm.dmFields = DM_PELSWIDTH | DM_PELSHEIGHT;
    
    ChangeDisplaySettings(&dm, CDS_FULLSCREEN);
    
    ShowWindow(hwnd, TRUE);
	UpdateWindow(hwnd);
    
    PIXELFORMATDESCRIPTOR pfd =
	{
		sizeof(PIXELFORMATDESCRIPTOR),
		1,
		PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER,    //Flags
		PFD_TYPE_RGBA,        // The kind of framebuffer. RGBA or palette.
		32,                   // Colordepth of the framebuffer.
		0, 0, 0, 0, 0, 0,
		0,
		0,
		0,
		0, 0, 0, 0,
		24,                   // Number of bits for the depthbuffer
		8,                    // Number of bits for the stencilbuffer
		0,                    // Number of Aux buffers in the framebuffer.
		PFD_MAIN_PLANE,
		0,
		0, 0, 0
	};

	HDC hdc = GetDC(hwnd);

	int  pf = ChoosePixelFormat(hdc, &pfd);
	SetPixelFormat(hdc, pf, &pfd);

	HGLRC glrc = wglCreateContext(hdc);
	wglMakeCurrent(hdc, glrc);
    
    glCreateShader = (PFNGLCREATESHADERPROC) wglGetProcAddress("glCreateShader");
	glCreateProgram = (PFNGLCREATEPROGRAMPROC) wglGetProcAddress("glCreateProgram");
	glShaderSource = (PFNGLSHADERSOURCEPROC) wglGetProcAddress("glShaderSource");
	glCompileShader = (PFNGLCOMPILESHADERPROC) wglGetProcAddress("glCompileShader");
	glAttachShader = (PFNGLATTACHSHADERPROC) wglGetProcAddress("glAttachShader");
	glLinkProgram = (PFNGLLINKPROGRAMPROC) wglGetProcAddress("glLinkProgram");
	glUseProgram = (PFNGLUSEPROGRAMPROC) wglGetProcAddress("glUseProgram");
    glGetUniformLocation = (PFNGLGETUNIFORMLOCATIONPROC) wglGetProcAddress("glGetUniformLocation");
    glUniform1f = (PFNGLUNIFORM1FPROC) wglGetProcAddress("glUniform1f");
    
    ShowCursor(FALSE);
    
    const char * gfx_source = "\
    #version 130\n\
    uniform float iTime;\
    const vec2 iResolution = vec2(1280.,720.);\
    void main()\
    {\
        vec2 uv = gl_FragCoord.xy/iResolution.xy;\
        vec3 col = 0.5 + 0.5*cos(iTime+uv.xyx+vec3(0,2,4));\
        gl_FragColor = vec4(col,1.0);\
    }";
    int gfx_size = strlen(gfx_source),
        gfx_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(gfx_handle, 1, (GLchar **)&gfx_source, &gfx_size);
    glCompileShader(gfx_handle);
    
    int gfx_program = glCreateProgram();
    glAttachShader(gfx_program, gfx_handle);
    glLinkProgram(gfx_program);
    
    glUseProgram(gfx_program);
    int gfx_iTime_location = glGetUniformLocation(gfx_program, "iTime");
    
    float t_now = 0.;
    
    while(1)
    {
        MSG msg = { 0 };
        while ( PeekMessageA( &msg, NULL, 0, 0, PM_REMOVE ) )
        {
            if ( msg.message == WM_QUIT || (msg.message = WM_KEYDOWN && msg.wParam == VK_ESCAPE ) )
            {
                ExitProcess(0);
                return 0;
            }
            
            TranslateMessage( &msg );
            DispatchMessageA( &msg );
        }
        
        glUniform1f(gfx_iTime_location, t_now);
        
        glViewport(0,0,1280,720);
        
        glBegin(GL_QUADS);
        glVertex3f(-1,-1,0);
        glVertex3f(-1,1,0);
        glVertex3f(1,1,0);
        glVertex3f(1,-1,0);
        glEnd();
        
        SwapBuffers(hdc);
    }
    
    return 0;
}
