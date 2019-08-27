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
 
#version 130

uniform float iTime;

vec2 uv,
    s,
    iResolution = vec2(1280.,720.);
vec3 c = vec3(1.,0.,-1.),
    o,
    dir,
    x,
    n,
    col,
    l;
float d = 0.,
    pi = acos(-1.);
int N = 550,
    i,
    j;

// Creative Commons Attribution-ShareAlike 4.0 International Public License
// Created by David Hoskins.
// See https://www.shadertoy.com/view/4djSRW
float hash12(vec2 p)
{
	vec3 p3  = fract(vec3(p.xyx) * .1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}

float lfnoise12(vec2 t)
{
    vec2 i = floor(t);
    t = fract(t);
    t = smoothstep(c.yy, c.xx, t);
    vec2 v1 = vec2(hash12(i), hash12(i+c.xy)), 
        v2 = vec2(hash12(i+c.yx), hash12(i+c.xx));
    v1 = c.zz+2.*mix(v1, v2, t.y);
    return mix(v1.x, v1.y, t.x);
}

float mfnoise12(vec2 x, float d, float b, float e)
{
    float n = 0.,
        a = 1., 
        nf = 0., 
        buf;
    for(float f = d; f<b; f *= 2.)
    {
        buf = lfnoise12(f*x);
        n += a*buf;
        a *= e;
        nf += 1.;
    }
    n *= (1.-e)/(1.-pow(e, nf));
    return n;
}

vec2 add(vec2 sda, vec2 sdb)
{
    return (sda.x<sdb.x)?sda:sdb;
}

vec2 scene(in vec3 y)
{
    y.y += .7*iTime;
    float da = mfnoise12(y.xy, 4.,450.,mix(.25,.5,smoothstep(.05,.35,y.z))),
    	db = da,
        r;
    
    da = .5 + .5*da;
    da *= smoothstep(.7,-.1,abs(y.x));
    da += .1*db;
    vec2 sdf = vec2(y.z-.6*da,1.),
        z = mod(y.xy, .4)-.2,
        zi = y.xy-z;
    r = hash12(zi);
    if(r > .9) 
    {
        sdf = add(sdf, vec2(abs(length(vec3(z,y.z-r+.9))-2.*(r-.9))-.001, 2.));
    }
    
    vec2 a = vec2(mod(y.x,.05)-.025, y.z);
    sdf = add(sdf, vec2(length(a)-.005, 3.));
    sdf = add(sdf, vec2(abs(length(a)-.025)-.001, 4.));
    return sdf;
}

vec3 normal(in vec3 y)
{
    float dx = 1.e-4,
        s = scene(y).x;
    return normalize(vec3(scene(y+dx*c.xyy).x-s, scene(y+dx*c.yxy).x-s, scene(y+dx*c.yyx).x-s));
}

vec3 colorize()
{
    float dln = dot(l,n),
        drev = abs(dot(reflect(normalize(x-l),n),dir));
    if(s.y == 1.)
    {
        return mix(.7*vec3(0.38,0.39,0.61),vec3(0.22,0.01,0.20),smoothstep(.0,.5,x.z))
            + .1*vec3(0.56,0.12,0.47)*dln
            + c.xxx*pow(drev,6.);
    }
    else if(s.y == 2.)
    {
        return .2*c.yyx
            + .1*c.xxx*dln
            + 8.*vec3(0.78,0.39,0.31)*pow(drev,2.);
    }
    else if(s.y == 3.)
    {
        return .1*vec3(1.00,0.87,0.74)
            + .1*vec3(1.00,0.87,0.74)*dln
            + 2.*vec3(1.00,0.87,0.74)*pow(drev,2.);
    }
    else if(s.y == 4.)
    {
        return c.xyy
            + c.xyy*dln
            + c.xyy*pow(drev,2.);
    }
}

void main()
{
    uv = (gl_FragCoord.xy-.5*iResolution.xy)/iResolution.yy;
    o = vec3(uv,1.);
    dir = c.yyz;
	float p = pi/4.;
    vec2 cs = vec2(cos(p), sin(p));
    mat2 R = mat2(cs.x,cs.y,-cs.y,cs.x);
    dir.yz = R * dir.yz;
    
    for(i=0; i<N; ++i)
    {
        x = o + d * dir;
        s = scene(x);
        if(s.x < 1.e-4) break;
        d += min(s.x,1.e-2);
        //d += s.x;
    }
    
    if(i<N)
    {
        n = normal(x);
        l = x + .1*n;
        col = colorize();
        
        if( s.y == 2. || s.y == 4.)
        {
            for(j = 0; j<2; ++j)
            {
                o = x;
                dir = refract(dir, n, .5);
                d = .1;

                for(i=0; i<N; ++i)
                {
                    x = o + d * dir;
                    s = scene(x);
                    if(s.x < 1.e-4) break;
                    d += min(s.x,5.e-2);
                }

                if(i<N)
                {
                    n = normal(x);
                    l = x + .1*n;
                    vec3 c1 = colorize();
                    col = mix(col, c1, .7);
                }
            }
        }
    }
    
    col *= col;
    
    gl_FragColor = vec4(clamp(col,0.,1.),1.0);
}
