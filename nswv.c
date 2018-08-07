/* 
 *  Eternal Darkness - 64k intro by QM^NR4/Team210
 * 
 *  Copyright (C) 2017  Alexander Kraus <nr4@z10.info>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <pthread.h>
#include <dlfcn.h>

#include <X11/X.h>
#include <X11/Xlib.h>
#include <X11/Xatom.h>
// #include <X11/Xutil.h>

#include <GL/gl.h>
#include <GL/glx.h>

#include <alsa/asoundlib.h>
#define PCM_DEVICE "default"

typedef void (*glGetShaderiv_t)(GLuint,  GLenum,  GLint *);
glGetShaderiv_t glGetShaderiv;
typedef void (*glGetShaderInfoLog_t)(GLuint,  GLsizei, GLsizei,  GLchar *);
glGetShaderInfoLog_t glGetShaderInfoLog;
typedef GLuint (*glCreateShader_t)(GLenum);
glCreateShader_t glCreateShader;
typedef GLuint (*glCreateProgram_t)();
glCreateProgram_t glCreateProgram;
typedef void (*glShaderSource_t)(GLuint, GLsizei, GLchar **, GLint *);
glShaderSource_t glShaderSource;
typedef void (*glCompileShader_t)(GLuint);
glCompileShader_t glCompileShader;
typedef void (*glAttachShader_t)(GLuint, GLuint);
glAttachShader_t glAttachShader;
typedef void (*glLinkProgram_t)(GLuint);
glLinkProgram_t glLinkProgram;
typedef void (*glUseProgram_t)(GLuint);
glUseProgram_t glUseProgram;
typedef GLint (*glGetUniformLocation_t)(GLuint, const GLchar *);
glGetUniformLocation_t glGetUniformLocation;
typedef void (*glUniform2f_t)(GLint, GLfloat, GLfloat);
glUniform2f_t glUniform2f;
typedef void (*glUniform1f_t)(GLint, GLfloat);
glUniform1f_t glUniform1f;
typedef void (*glGenFramebuffers_t)(GLsizei, GLuint*);
glGenFramebuffers_t glGenFramebuffers;
typedef void (*glBindFramebuffer_t)(GLenum, GLuint);
glBindFramebuffer_t glBindFramebuffer;
typedef void (*glFramebufferTexture2D_t)(GLenum, GLenum, GLenum, GLuint, GLint);
glFramebufferTexture2D_t glFramebufferTexture2D;
typedef void (*glBlitFramebuffer_t)(GLint, GLint, GLint, GLint, GLint, GLint, GLint, GLint, GLbitfield, GLenum);
glBlitFramebuffer_t glBlitFramebuffer;
typedef void (*glNamedRenderbufferStorage_t) (GLuint, GLenum, GLsizei, GLsizei);
glNamedRenderbufferStorage_t glNamedRenderbufferStorage;

Display                 *display;
Window                  root;
GLint                   att[] = { GLX_RGBA, GLX_DEPTH_SIZE, 24, GLX_DOUBLEBUFFER, None };
XVisualInfo             *vi;
Colormap                cmap;
XSetWindowAttributes    swa;
Window                  win;
GLXContext              glc;
XWindowAttributes       gwa;
XEvent                  xevent;

snd_pcm_t *snd_handle;
snd_pcm_hw_params_t *params;
snd_pcm_uframes_t frames;

struct timeval t;

unsigned int sample_rate = 44100, channels = 2;
double duration1 = 112./77.5*60./2., duration2 = 100.;
float *music1, *smusic1, *music2, *smusic2;
int music1_size, music2_size;

int scene_index = 0;
double newtime = .5;

double volume = .1;

//FIXME DEBUG ONLY
void debug(int shader_handle)
{
    int compile_status = 0;
    glGetShaderiv(shader_handle, GL_COMPILE_STATUS, &compile_status);
    if(compile_status != GL_TRUE)
    {
        printf("FAILED.\n");
        int len = 12;
        glGetShaderiv(shader_handle, GL_INFO_LOG_LENGTH, &len);
        printf("log length: %d\n", len);
        GLchar CompileLog[len];
        glGetShaderInfoLog(shader_handle, len, NULL, CompileLog);
        printf("error: %s\n", CompileLog);
    }
}
//FIXME END DEBUG ONLY

// void print32(float *f)
// {
//     int d = *(int*)f;
//     for(int i=0; i<32; ++i)
//     {
//         int mask = 1<<i;
// //         printf("%d\n", mask);
//         int bit = d & mask;
//         if(bit != 0) printf("1");
//         else printf("0");
//         if(i%8 == 7) printf(" ");
//     }
//     printf("\n");
// }

void separate(float *in, float *o1, float *o2)
{
    unsigned int d = *(int*)in, low_bytes = 65535, high_bytes = 65535<<16;
    float *hif = (float*)&low_bytes, *lof = (float*)&high_bytes;
//     print32(hif);
//     printf("low mask: ");
//     print32(lof);
    
    unsigned int d1 = (d & low_bytes) /*<< 16*/;
    float *d1f = (float*)&d1;
    
//     printf("low bytes: ");
//     print32(d1f);
    
    unsigned int d2 = (d & high_bytes) >> 16;
    float *d2f = (float*)&d2;
    
//     printf("high bytes: ");
//     print32(d2f);
    
    *o1 = d1;
    *o2 = d2;
    *o1 = 2.f**o1/(float)(low_bytes+1)-1.f;
    *o2 = 2.f**o2/(float)(low_bytes+1)-1.f;
//     printf("extracted: %le %le\n", 2.f**o1/(float)(low_bytes+1)-1.f, 2.f**o2/(float)(low_bytes+1)-1.f);
//     printf("
}

void playMusic()
{
    float *d = (float*)malloc(10*sample_rate*sizeof(float));
    for(int i=0; i<10*sample_rate; ++i)
        d[i] = 0.;
    snd_pcm_writei(snd_handle, d, 5*sample_rate);
//     snd_pcm_drain(snd_handle);
//     printf("%le: finished silence\n", newtime);
    free(d);
//     while(newtime <5.);
    snd_pcm_writei(snd_handle, smusic1, 2.*(music1_size+music2_size));
//     snd_pcm_drain(snd_handle);
//     printf("%le: finished no1\n", newtime);
//     while(newtime< 10.);

//     while(snd_pcm_drain(snd_handle)<0);
// //     snd_pcm_writei(snd_handle, smusic2, 2.*music2_size);
//     printf("%le: finished no2\n", newtime);

}

void letexit()
{
    while(1)
    {
        while(XPending(display))
        {
            XNextEvent(display, &xevent);
        }
        if(xevent.type == KeyPress) 
        {
            exit(0);
        }
    }
}

int main(int argc, char **args)
{
    void *gl = (void*)dlopen("libGL.so", RTLD_LAZY | RTLD_GLOBAL);
    
    glGetShaderiv = (glGetShaderiv_t) dlsym(gl, "glGetShaderiv");
    glGetShaderInfoLog = (glGetShaderInfoLog_t) dlsym(gl, "glGetShaderInfoLog");
    glCreateShader = (glCreateShader_t) dlsym(gl, "glCreateShader");
    glCreateProgram = (glCreateProgram_t) dlsym(gl, "glCreateProgram");
    glShaderSource = (glShaderSource_t) dlsym(gl, "glShaderSource");
    glCompileShader = (glCompileShader_t) dlsym(gl, "glCompileShader");
    glAttachShader = (glAttachShader_t) dlsym(gl, "glAttachShader");
    glLinkProgram = (glLinkProgram_t) dlsym(gl, "glLinkProgram");
    glUseProgram = (glUseProgram_t) dlsym(gl, "glUseProgram");
    glGetUniformLocation = (glGetUniformLocation_t) dlsym(gl, "glGetUniformLocation");
    glUniform2f = (glUniform2f_t) dlsym(gl, "glUniform2f");
    glUniform1f = (glUniform1f_t) dlsym(gl, "glUniform1f");
    glGenFramebuffers = (glGenFramebuffers_t) dlsym(gl, "glGenFramebuffers");
    glBindFramebuffer = (glBindFramebuffer_t) dlsym(gl, "glBindFramebuffer");
    glFramebufferTexture2D = (glFramebufferTexture2D_t) dlsym(gl, "glFramebufferTexture2D");
    glBlitFramebuffer = (glBlitFramebuffer_t) dlsym(gl, "glBlitFramebuffer");
    glNamedRenderbufferStorage = (glNamedRenderbufferStorage_t) dlsym(gl, "glNamedRenderbufferStorage");
    
    if(gl == 0) 
    {
        printf("COULD NOT DLOPEN LIBGL. EXPLOSION NOW.\n");
        printf("%s\n", dlerror());
    }
    
    XInitThreads();
    
    //create window and init opengl
    display = XOpenDisplay(NULL);
    root = DefaultRootWindow(display);
    vi = glXChooseVisual(display, 0, att);
    cmap = XCreateColormap(display, root, vi->visual, AllocNone);
    swa.colormap = cmap;
    int w=1920, h=1080;
    if(argc >= 6)
    {
        sscanf(args[5], "%le", &volume);
    }
    printf("Playing with volume: %le\n", volume);
    if(argc >= 3)
    {
        sscanf(args[1], "%d", &w);
        sscanf(args[2], "%d", &h);
    }
//     printf(
    swa.event_mask = ExposureMask | KeyPressMask;
    win = XCreateWindow(display, root, 0, 0, w, h, 0, vi->depth, InputOutput, vi->visual, CWColormap | CWEventMask, &swa);
 
    Atom atoms[2] = { XInternAtom(display, "_NET_WM_STATE_FULLSCREEN", True), None };
    XChangeProperty(
        display, 
        win, 
        XInternAtom(display, "_NET_WM_STATE", True),
                    XA_ATOM, 32, PropModeReplace,(unsigned char*) atoms, 1
    );
    XMapWindow(display, win);
    glc = glXCreateContext(display, vi, NULL, GL_TRUE);
    glXMakeCurrent(display, win, glc);
    
//     glewInit();
    
    pthread_t *exithread;
    pthread_create(&exithread, NULL, (void*)&letexit, 0);
    
    int maxUniformVectors;
    glGetIntegerv(GL_MAX_FRAGMENT_UNIFORM_VECTORS, &maxUniformVectors);
    printf("maxuniformvectors: %d\n", maxUniformVectors);
    
    //display loading bar
#include "load210.h"
    printf("LOADING BAR\n");
    int lb_size = strlen(load210_frag), lb_handle = glCreateShader(GL_FRAGMENT_SHADER), lb_program = glCreateProgram();
    glShaderSource(lb_handle, 1, (const char *const *)&load210_frag, &lb_size);
    glCompileShader(lb_handle);
    debug(lb_handle);
    glAttachShader(lb_program, lb_handle);
    glLinkProgram(lb_program);
    glUseProgram(lb_program);
    int lb_progress_location = glGetUniformLocation(lb_program, VAR_PROGRESS),
        lb_time_location = glGetUniformLocation(lb_program, VAR_ITIME),
        lb_resolution_location = glGetUniformLocation(lb_program, VAR_IRESOLUTION);
//     XGetWindowAttributes(display, win, &gwa);
    glUniform2f(lb_resolution_location, w, h);
    
    //open pcm device and setup asoundlib
    int block_size = 512*512;
    int nblocks1 = (int)ceil((double)(sample_rate*channels)*duration1/(double)block_size);
    int nblocks2 = (int)ceil((double)(sample_rate*channels)*duration2/(double)block_size);
    music1_size = (nblocks1+nblocks2)*block_size; 
    music1 = (float*)malloc(music1_size*sizeof(float)),
    smusic1 = (float*)malloc(2*music1_size*sizeof(float));
//     music2_size = nblocks2*block_size; 
//     music2 = (float*)malloc(music2_size*sizeof(float)),
//     smusic2 = (float*)malloc(2*music2_size*sizeof(float));
    
    double newtime = .5;
        
        glUseProgram(lb_program);
        
        glUniform1f(lb_time_location, newtime);
        glUniform1f(lb_progress_location, 1.);
        glUniform2f(lb_resolution_location, w, h);
        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glViewport(0,0,w,h);
        glBegin(GL_QUADS);
        glVertex3f(-1,-1,0);
        glVertex3f(-1,1,0);
        glVertex3f(1,1,0);
        glVertex3f(1,-1,0);
        glEnd();
        
        glXSwapBuffers(display, win);
    
    glUseProgram(0);
    
    
    
    snd_pcm_open(&snd_handle, PCM_DEVICE, SND_PCM_STREAM_PLAYBACK, 0);
    snd_pcm_set_params(snd_handle, SND_PCM_FORMAT_FLOAT, SND_PCM_ACCESS_RW_INTERLEAVED, channels, sample_rate, 0, (music1_size+music2_size)*channels);
    
    //generate textures for music1 rendering
    unsigned int snd_framebuffer = 0;
    glGenFramebuffers(1, &snd_framebuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, snd_framebuffer);
    
    unsigned int snd_texture;
    glGenTextures(1, &snd_texture);
    glBindTexture(GL_TEXTURE_2D, snd_texture);
    glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA, 512, 512, 0,GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, snd_texture, 0);

#undef VAR_ITIME
#undef VAR_IRESOLUTION
#include "skyscrapers.h"
    printf("SKYSCRAPERS\n");
    int ss_size = strlen(skyscrapers_frag), ss_handle = glCreateShader(GL_FRAGMENT_SHADER), ss_program = glCreateProgram();
    glShaderSource(ss_handle, 1, (const char *const *)&skyscrapers_frag, &ss_size);
    glCompileShader(ss_handle);
    debug(ss_handle);
    glAttachShader(ss_program, ss_handle);
    glLinkProgram(ss_program);
    glUseProgram(ss_program);
    int ss_time_location = glGetUniformLocation(ss_program, VAR_ITIME),
        ss_resolution_location = glGetUniformLocation(ss_program, VAR_IRESOLUTION),
        ss_intensity_location = glGetUniformLocation(ss_program, VAR_INTENSITY);
//     XGetWindowAttributes(display, win, &gwa);
//     glUniform2f(lb_resolution_location, w, h);
    
    glUseProgram(0);
    
#undef VAR_ITIME
#undef VAR_IRESOLUTION
#include "chipscene.h"
    printf("CHIPS\n");
    int ch_size = strlen(chipscene_frag), ch_handle = glCreateShader(GL_FRAGMENT_SHADER), ch_program = glCreateProgram();
    glShaderSource(ch_handle, 1, (const char *const *)&chipscene_frag, &ch_size);
    glCompileShader(ch_handle);
    debug(ch_handle);
    glAttachShader(ch_program, ch_handle);
    glLinkProgram(ch_program);
    glUseProgram(ch_program);
    int ch_time_location = glGetUniformLocation(ch_program, VAR_ITIME),
        ch_resolution_location = glGetUniformLocation(ch_program, VAR_IRESOLUTION),
        ch_intensity_location = glGetUniformLocation(ch_program, VAR_INTENSITY);
//     XGetWindowAttributes(display, win, &gwa);
    glUniform2f(ch_resolution_location, w, h);
    
    glUseProgram(0);

    
    //generate audio with audio shader
#include "sound.h"
    printf("SOUND\n");
    
    int as_size = strlen(sound_frag), as_handle = glCreateShader(GL_FRAGMENT_SHADER), as_program = glCreateProgram();
    glShaderSource(as_handle, 1, (const char *const *)&sound_frag, &as_size);
    glCompileShader(as_handle);
    debug(as_handle);
    glAttachShader(as_program, as_handle);
    glLinkProgram(as_program);
    glUseProgram(as_program);
    int as_sample_rate_location = glGetUniformLocation(as_program, VAR_ISAMPLERATE),
        as_block_offset_location = glGetUniformLocation(as_program, VAR_IBLOCKOFFSET),
        as_ivolume_location = glGetUniformLocation(as_program, VAR_IVOLUME);
    glUseProgram(0);

#include "sound2.h"
    printf("SOUND2\n");
    
    int as2_size = strlen(sound2_frag), as2_handle = glCreateShader(GL_FRAGMENT_SHADER), as2_program = glCreateProgram();
    glShaderSource(as2_handle, 1, (const char *const *)&sound2_frag, &as2_size);
    glCompileShader(as2_handle);
    debug(as2_handle);
    glAttachShader(as2_program, as2_handle);
    glLinkProgram(as2_program);
    glUseProgram(as2_program);
    int as2_sample_rate_location = glGetUniformLocation(as2_program, VAR_ISAMPLERATE),
        as2_block_offset_location = glGetUniformLocation(as2_program, VAR_IBLOCKOFFSET),
        as2_ivolume_location = glGetUniformLocation(as2_program, VAR_IVOLUME);
    glUseProgram(0);
    
    
    
    struct timeval tv;
    gettimeofday(&tv, NULL);
    unsigned long long t0 = (unsigned long long)(tv.tv_sec) * 1000 + (unsigned long long)(tv.tv_usec) / 1000;
//     double newtime;
    
    int x_file_descriptor = ConnectionNumber(display);
    fd_set x_file_descriptor_set;
    for(int i=0; i<nblocks1; ++i)
    {
        double tstart = (double)(i*block_size)/(double)sample_rate;
        
        glBindFramebuffer(GL_FRAMEBUFFER, snd_framebuffer);
        glViewport(0,0,512,512);
        
        glUseProgram(as_program);
        glUniform1f(as_ivolume_location, (float)volume);
        glUniform1f(as_sample_rate_location, (float)sample_rate);
        glUniform1f(as_block_offset_location, (float)tstart);
        
//         glDrawBuffers(1, DrawBuffers);

        glBegin(GL_QUADS);
            glVertex3f(-1,-1,0);
            glVertex3f(-1,1,0);
            glVertex3f(1,1,0);
            glVertex3f(1,-1,0);
        glEnd();
        
        glReadPixels(0, 0, 512, 512, GL_RGBA, GL_UNSIGNED_BYTE, music1+i*block_size);
        
        glUseProgram(0);
        
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        
        for(int j=0; j<block_size; ++j)
        {
            separate(music1+i*block_size+j, smusic1+2*i*block_size+2*j, smusic1+2*i*block_size+2*j+1);
        }
        
        newtime = .5;
        
        glUseProgram(lb_program);
        
        glUniform1f(lb_time_location, newtime);
        glUniform1f(lb_progress_location, 1.-(float)i/(float)(nblocks1+nblocks2));
        glUniform2f(lb_resolution_location, w, h);
        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glViewport(0,0,w,h);
        glBegin(GL_QUADS);
        glVertex3f(-1,-1,0);
        glVertex3f(-1,1,0);
        glVertex3f(1,1,0);
        glVertex3f(1,-1,0);
        glEnd();
        
        glXSwapBuffers(display, win);
        
        glUseProgram(0);
        
//         printf("Done Block %d: tstart = %le\n", i, tstart);
        
    }

    for(int i=nblocks1; i<nblocks1+nblocks2; ++i)
    {
//         printf("started lloop %d\n", i-nblocks1);
        double tstart = (double)((i-nblocks1)*block_size)/(double)sample_rate;
//         printf("tstart = %le\n", tstart);
        glBindFramebuffer(GL_FRAMEBUFFER, snd_framebuffer);
        glViewport(0,0,512,512);
        
        glUseProgram(as2_program);
        
        glUniform1f(as2_ivolume_location, (float)volume);
        glUniform1f(as2_sample_rate_location, (float)sample_rate);
        glUniform1f(as2_block_offset_location, (float)tstart);
        
//         glDrawBuffers(1, DrawBuffers);

        glBegin(GL_QUADS);
            glVertex3f(-1,-1,0);
            glVertex3f(-1,1,0);
            glVertex3f(1,1,0);
            glVertex3f(1,-1,0);
        glEnd();
        
//         printf("trying to read %d from graphics\n", i-nblocks1);
        glReadPixels(0, 0, 512, 512, GL_RGBA, GL_UNSIGNED_BYTE, music1+i*block_size);
//         printf("read %d from graphics\n", i-nblocks1);
        glUseProgram(0);
        
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        
        for(int j=0; j<block_size; ++j)
        {
            separate(music1+nblocks1*block_size+(i-nblocks1)*block_size+j, smusic1+2*nblocks1*block_size+2*(i-nblocks1)*block_size+2*j, smusic1+2*nblocks1*block_size+2*(i-nblocks1)*block_size+2*j+1);
        }
//         printf("separated %d\n", i-nblocks1);
        newtime = .5;
        
        glUseProgram(lb_program);
        
        glUniform1f(lb_time_location, newtime);
        glUniform1f(lb_progress_location, 1.-(float)i/(float)(nblocks1+nblocks2));
        glUniform2f(lb_resolution_location, w, h);
        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glViewport(0,0,w,h);
        glBegin(GL_QUADS);
        glVertex3f(-1,-1,0);
        glVertex3f(-1,1,0);
        glVertex3f(1,1,0);
        glVertex3f(1,-1,0);
        glEnd();
        
        glXSwapBuffers(display, win);
        
        glUseProgram(0);
        
//         printf("Done Block %d: tstart = %le\n", i, tstart);
        
    }

//     unsigned int sky_framebuffer = 0;
//     glGenFramebuffers(1, &sky_framebuffer);
//     glBindFramebuffer(GL_FRAMEBUFFER, sky_framebuffer);
    
//     unsigned int sky_texture;
//     glGenTextures(1, &sky_texture);
//     glBindTexture(GL_TEXTURE_2D, sky_texture);
//     glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA, 640, 480, 0,GL_RGBA, GL_UNSIGNED_BYTE, 0);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

//     glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, sky_texture, 0);
//     for(int i=0; i<2.*music1_size; ++i)
//         printf("%le %le\n", (float)i/(float)sample_rate, smusic1[i]);
    


//     GLenum DrawBuffers[1] = {GL_COLOR_ATTACHMENT0};
    
    pthread_t music1_thread;
    pthread_create(&music1_thread, NULL, (void*)&playMusic, 0);
    
//     glDeleteFramebuffers(1, &snd_framebuffer);
    
    int usedxres = 860, usedyres=540;
    if(argc >= 5)
    {
        sscanf(args[3], "%d", &usedxres);
        sscanf(args[4], "%d", &usedyres);
    }
    printf("using xres: %d %d\n", usedxres, usedyres);
    printf("using res: %d %d\n", w, h);
    glUseProgram(ss_program);
    glUniform2f(ss_resolution_location, usedxres, usedyres);
//     GLuint renderbuffer;
//     glCreateRenderbuffers(1, &renderbuffer);
    
//     glFramebufferRenderbuffer(snd_framebuffer, GL_COLOR_ATTACHMENT0, 
    
//     glNamedRenderbufferStorage(snd_framebuffer, GL_RGBA, w, h);
    glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA, w, h, 0,GL_RGBA, GL_UNSIGNED_BYTE, 0);
    
    double dt = newtime;
    //=============== RENDER LOOP ================
    t0 = (unsigned long long)(tv.tv_sec) * 1000 + (unsigned long long)(tv.tv_usec) / 1000;
    while(1) 
    {
        
//         if(xevent.type == KeyPress) 
//         {
//             return 0;
//         }
        
        FD_ZERO(&x_file_descriptor_set);
        FD_SET(x_file_descriptor, &x_file_descriptor_set);
        
        t.tv_usec = 1.e6/60.;
        t.tv_sec = 0;
        
        int num_ready_fds = select(x_file_descriptor + 1, &x_file_descriptor_set, NULL, NULL, &t);
        if (num_ready_fds == 0)
        {
            
            //TODO: add loading mechanism here.
            gettimeofday(&tv, NULL);
            unsigned long long newtime_l = (unsigned long long)(tv.tv_sec) * 1000 + (unsigned long long)(tv.tv_usec) / 1000 - t0;
            double newtime = (double)newtime_l/1000.;
            
//             printf("%le\n", newtime);
            double tmusic1 = newtime-10.25;
            double scale = 0.;
            
            if(tmusic1 >= 2.*duration1+1.5) scene_index = 1.;
            if(tmusic1 > 1.)
            {
                int l = 5500;
                for(int i=-l; i<l; ++i)
                {
//                     printf("%d\n", (int)tmusic1*sample_rate*channels+i);
                    scale += fabs(smusic1[(int)(tmusic1*(double)sample_rate*(double)channels)+i]);
                }
                scale = scale /2./(double)l;
//                 printf("%le\n", scale);
                
            }
            if(tmusic1 >= 2*duration1+2*duration2+.3)exit(0);
            if(scene_index == 0)
            {
                
                glUseProgram(ss_program);
                
                glUniform1f(ss_intensity_location, scale);
                glUniform1f(ss_time_location, newtime);
                
                glBindFramebuffer(GL_DRAW_FRAMEBUFFER, snd_framebuffer);
                glViewport(0,0,usedxres,usedyres);
                
                glClearColor(0.,0.,0.,1.);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                
                glBegin(GL_QUADS);
                glVertex3f(-1,-1,0);
                glVertex3f(-1,1,0);
                glVertex3f(1,1,0);
                glVertex3f(1,-1,0);
                glEnd();
    //             glViewport(0,0,w,h);
                glBindFramebuffer(GL_READ_FRAMEBUFFER, snd_framebuffer);
                glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    //             glViewport(0,0, w, h);
                glBlitFramebuffer(0,0,usedxres,usedyres,0,0,w,h, GL_COLOR_BUFFER_BIT, GL_LINEAR);
            }
            else if(scene_index == 1)
            {
                glUseProgram(ch_program);
                glUniform1f(ch_intensity_location, scale);
                
                glViewport(0,0,w,h);
                glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
                glUniform1f(ch_time_location, newtime);
                glClearColor(0.,0.,0.,1.);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                
                glBegin(GL_QUADS);
                glVertex3f(-1,-1,0);
                glVertex3f(-1,1,0);
                glVertex3f(1,1,0);
                glVertex3f(1,-1,0);
                glEnd();
            }
            
            glXSwapBuffers(display, win);
            glUseProgram(0);
        }
    }
    
    dlclose("libGL.a");
    return 0;
}
