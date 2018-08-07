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

#version 130

#define AA 2

const vec2 c = vec2(1.,0.);
const float pi = acos(-1.);

uniform float progress, iTime;
uniform vec2 iResolution;
out vec4 fragCoord;

float rand(vec2 a0) 
{
    return fract(sin(dot(a0.xy ,vec2(12.9898,78.233)))*43758.5453);
}

float blend(float a) 
{
    return ((6.*a-15.)*a+10.)*a*a*a;
}

float g(float a, float a0) 
{
    return (2.*rand(vec2(a,a0))-1.);
}

float perlin1d(float a, float seed) 
{
    return mix(g(floor(a),seed)*fract(a),dot(vec2(g(floor(a)+1.,seed)),vec2(fract(a),-1.)),blend(fract(a)));
}

float mfperlin1d(float x, float seed, float fmin, float fmax, float phi)
{
    float sum = 0.;
    float a = 1.;
    
    for(float f = fmin; f<fmax; f = f*2.)
    {
        sum = a*perlin1d(f*x, seed) + sum;
        a = a*phi;
    }
    
    return sum;
}

#define circle(a,a0) (length(a)-(a0))
#define box(a,a0,a1) (length(max(abs(a-(a0))-(a0),0.))-(a1))

#define add(a,a0) mix(a0,a,step((a).x,(a0).x))
#define sec(a,a0) mix(a,-(a0,step((a).x,(a0).x))

vec2 scene(vec2 x)
{
    float phi = -.75*pi;
    mat2 m = mat2(sin(phi), cos(phi), -cos(phi), sin(phi));
    
    vec2 sdf = c.xy;
    
    //210
    sdf = add(sdf, vec2(circle(x-vec2(.4,.4), .2),1.));
    sdf = add(sdf, vec2(circle(x-vec2(.7,.4), .2),1.));
	sdf = add(sdf, vec2(box(x-vec2(.5,.2), vec2(.05,.2), 0.), 1.));
    
    //T
    sdf = add(sdf, vec2(box(x-vec2(.45, .05), vec2(.01, .05), 0.), 2.));
    sdf = add(sdf, vec2(box(x-vec2(.42, .15), vec2(.04, .015), 0.), 2.));
    
    x.x -= .01;
    
    //E
    sdf = add(sdf, vec2(box(x-vec2(.5, .05), vec2(.01, .065), 0.), 2.));
    sdf = add(sdf, vec2(box(x-vec2(.5, .15), vec2(.03, .015), 0.), 2.));
    sdf = add(sdf, vec2(box(x-vec2(.5, .1), vec2(.03, .015), 0.), 2.));
    sdf = add(sdf, vec2(box(x-vec2(.5, .05), vec2(.03, .015), 0.), 2.));
    
    x.x -= .02;
    
    //A
    sdf = add(sdf, vec2(box(x-vec2(.55, .05), vec2(.01, .065), 0.), 2.));
    sdf = add(sdf, vec2(box(x-vec2(.6, .05), vec2(.01, .065), 0.), 2.));
    sdf = add(sdf, vec2(box(x-vec2(.55, .15), vec2(.03, .015), 0.), 2.));
    sdf = add(sdf, vec2(box(x-vec2(.55, .1), vec2(.03, .015), 0.), 2.));
    
    x.x -= .03;
    
    //M
    sdf = add(sdf, vec2(box(x-vec2(.6, .05), vec2(.01, .065), 0.), 2.));
    sdf = add(sdf, vec2(box(x-vec2(.63, .1), vec2(.01, .04), 0.), 2.));
    sdf = add(sdf, vec2(box(x-vec2(.66, .05), vec2(.01, .065), 0.), 2.));
    sdf = add(sdf, vec2(box(x-vec2(.6, .15), vec2(.04, .015), 0.), 2.));
    
    return sdf;
}

float bubble(vec2 x)
{
    x.y-= + mix(-.9,.9,progress);
    vec2 y = mod(x, .4*c.xx)-.2*c.xx;
    vec2 index = x-y;
        
    vec2 center = c.yy-.05+.05*rand(index);
  	float r = length(y-center), width = .0015, rc = .07*rand(index);
    float scale = 1./width/sqrt(2.*pi)*exp(-.5*pow((r-rc)/width, 2.));
    return scale;
}

vec3 colorize(vec2 x)
{
    vec3 color = c.yyy;
    
    vec2 xb = x;
    
    x += vec2(2.8, 1.2);
    x.xy*=.3;
    x.y-=.1*cos(2.4*pi*x.x-pi-3.*iTime);
    x.x-=.5*x.y;
    
    vec2 sdf = scene(x);
    float phi = atan(xb.y/xb.x);
    vec3 bg = mix(vec3(75.,52.,232.)/255., vec3(75.,138.,232.)/255., round(1.*(.5+.5*sin(6.*phi+3.e-1*iTime)))/1.);
    
    if(sdf.x > 0. || x.x <.4) 
        return bg;
    
    if(sdf.y < 0.) return c.yyy;
    
    if(sdf.y == 1.)
    {
        
        if(abs(sdf.x) < .1) return c.yyy;
		//if(x.x < .5) return bg;
        
        color = 0.5 + 0.5*cos(iTime+xb.xxy+vec3(0,2,4));
        
        float scale = bubble(xb);
        for(float i=0.; i<5.; i+=1.)
	        scale = max(scale, bubble(xb+2.*vec2(rand(vec2(i,i+1.)),rand(vec2(i, i-1.)))));

        color = mix(color, .4*c.xxx, clamp(scale, 0.,1.));
        
        if(xb.y < mix(-.9,.9,progress) + .18*mfperlin1d(xb.x-1.e-1*(4.*iTime), 0., 8.e0, 1.e3, .2) + .05*sin(10.*iTime-10.*xb.x))
	        color = mix(color, .4*c.xxx, clamp(scale, 0.,1.));
        else 
            color = mix(c.xxx, .5*c.xxx, clamp(scale, 0.,1.));
    }
    
    if(sdf.y == 2.)
    {
        return c.yyy;
    }
    
    return color;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec3 col = c.yyy;
    
    fragCoord *= 2.;
    fragCoord -= .5*iResolution;

    vec2 p;
    
#if AA!=1
    for(int i=0; i<AA; ++i)
    {
    	for(int j=0; j<AA; ++j)
        {
            vec2 o = vec2(float(i),float(j)) / float(AA) - 0.5;
        	p = (-iResolution.xy + 2.*(fragCoord+o))/iResolution.y;
#else 
            p = (-iResolution.xy + 2.*fragCoord)/iResolution.y;
#endif
            col += colorize(p);
            
#if AA!=1
        }
    }
    col/=float(AA*AA);
#else
#endif

    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
