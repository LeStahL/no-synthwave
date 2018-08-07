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

uniform float iTime;
uniform vec2 iResolution;
uniform float intensity;

#define AA 1

#define T .5

const vec2 c = vec2(1.,0.);
const float pi = 3.14159;
float t;

const float dxy = 1.;

float textw = .5;
float texth = 1.;

mat2 inverse(mat2 m)
{
    return mat2(m[1][1], -m[1][0], -m[0][1], m[0][0])/(m[0][0]*m[1][1]-m[1][0]*m[0][1]);
}

#define rand(a0) fract(sin(dot(a0.xy ,vec2(12.9898,78.233)))*43758.5453)
#define blend(a) ((6.*a-15.)*a+10.)*a*a*a
#define interpolate(d,w00,w10,w01,w11) mix(mix(w00,w10,blend(d.x)),mix(w01,w11,blend(d.x)),blend(d.y))

vec2 g2d(vec2 x, float seed)
{
    return vec2(-1.)+2.*vec2(rand(x+vec2(seed+2., seed+1.)), rand(x+vec2(seed+3.,seed+4.)));
}

float perlin2d(vec2 x, float seed)
{
    return interpolate(fract(x),
                       dot(g2d(floor(x), seed), fract(x)), 
                       dot(g2d(floor(x).xy+c.xy, seed), fract(x).xy-c.xy), 
                       dot(g2d(floor(x).xy+c.yx, seed), fract(x).xy-c.yx), 
                       dot(g2d(floor(x).xy+c.xx, seed), fract(x).xy-c.xx));
}

float mfperlin2d(vec2 x, float seed, float fmin, float fmax, float phi)
{
    float sum = 0.;
    float a = 1.;
    
    for(float f = fmin; f<fmax; f = f*2.)
    {
        sum = a*perlin2d(f*x, seed) + sum;
        a = a*phi;
    }
    
    return sum;
}

#define g(a,a0) (2.*rand(vec2(a,a0))-1.)
#define perlin1d(a,seed) mix(g(floor(a),seed)*fract(a),dot(vec2(g(floor(a)+1.,seed)),vec2(fract(a),-1.)),blend(fract(a)))
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

#define Rx(x) mat3(c.xyyy, cos(x), sin(x), 0., -sin(x), cos(x))
#define Ry(x) mat3(cos(x), 0., -sin(x), c.yxy, sin(x), 0., cos(x))
#define Rz(x) mat3(cos(x), sin(x), 0., -sin(x), cos(x), c.yyyx)
#define rotate(a,a0) Rz(-(a0).z)*Ry(-(a0).y)*Rx(-(a0).x)*(a)

mat3 rot(vec3 p)
{
    vec3 cp = cos(p), sp = sin(p);
    mat3 m = mat3(cp.y*cp.x, cp.x*sp.z+cp.z*sp.x*sp.y, sp.x*sp.z-cp.x*cp.z*sp.y, 
           -cp.y*sp.z, cp.x*cp.z-sp.x*sp.y*sp.z, cp.z*sp.x+cp.x*sp.y*sp.z, 
           sp.y, -cp.y*sp.x, cp.x*cp.y);
    return m;
}

#define sphere(a,a0) length(a)-(a0)
#define box(a,a0) length(max(abs(a)-(a0),0.))
float cylinder(vec3 x, vec2 h)
{
  vec2 d = abs(vec2(length(x.xz),x.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0))-.02;
}
float cylinder2(vec3 x, vec2 h)
{
  vec2 d = abs(vec2(length(x.xz),x.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}
#define sphere2(a,a0,a1) vec2(sphere(a,a0),a1)
#define box2(a,a0,a1) vec2(box(a,a0),a1)
#define plane2(a,a0,a1) vec2((a).z-(a0),a1)

#define repeat(a,a0) mod(a,a0)-.5*a0
#define add(a,a0) mix(a0,a,step((a).x,(a0).x))
#define sub(a,a0) mix(a0,a,step((a).x,-(a0).x))

// exponential smooth min (k = 32);
float sminexp( float a, float b, float k )
{
    float res = exp( -k*a ) + exp( -k*b );
    return -log( res )/k;
}

// polynomial smooth min (k = 0.1);
float sminpol( float a, float b, float k )
{
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

// power smooth min (k = 8);
float sminpow( float a, float b, float k )
{
    a = pow( a, k ); b = pow( b, k );
    return pow( (a*b)/(a+b), 1.0/k );
}

vec2 addexp(vec2 a, vec2 b)
{
    vec2 ret = add(a,b);
    ret.x = sminexp(a.x, b.x, 12.);
    return ret;
}

vec2 addpol(vec2 a, vec2 b)
{
    vec2 ret = add(a,b);
    ret.x = sminpol(a.x, b.x, .1);
    return ret;
}

vec2 addpow(vec2 a, vec2 b)
{
    vec2 ret = add(a,b);
    ret.x = sminpow(a.x, b.x, 8.);
    return ret;
}

vec3 bendx(vec3 x, float angle)
{
    float c = cos(angle*x.y);
    float s = sin(angle*x.y);
    mat2  m = mat2(c,-s,s,c);
    x = vec3(m*x.xy,x.z);
    return x;
}

vec2 chip(vec3 x, float color)
{
    vec2 sdf = c.xy;
    
    float delta = .04;
    
    sdf = add(sdf, vec2(box(x+1.*c.yyx, vec3(.2,.2,.1)), color));
    
    vec3 y = vec3(mod(x.xy,delta)-.5*delta, x.z+1.);
    
    vec2 index = (x-y).xy;
    if(length(index) < .25)
	    sdf = add(sdf, vec2(box(y, vec3(.2*delta,.2*delta,.1)), 3.));
    
    float guard = -box(y,delta*c.xxy+.2*c.yyx);
    guard = abs(guard)+delta*.2;
    sdf.x = min(sdf.x, guard);
    
    return sdf;
}

vec2 zap(vec3 x)
{
    vec2 sdf = c.xy;
    
    x.x += .8*mfperlin1d(x.z, iTime, 1., 1.e2, .4);
    x.y += .8*mfperlin1d(x.z, iTime, 1., 1.e2, .4);

    sdf = add(sdf, vec2(cylinder2(x.yzx, vec2(5.e-3, 100.)), 5.));
    
    return sdf;
}

vec2 scene(vec3 x)
{
    vec2 sdf = c.xy;
    //x = rot(vec3(1.,2.,3.)*1.e-1*iTime)*x;
    //x = bendx(x-.5*c.yxy, .3*sin(iTime));
    
    x += 1.*c.yxy+iTime*c.yxy;
    
    sdf = add(sdf, plane2(x,-1., 15.));
    
    vec3 y = vec3(mod(x.xy, dxy)-.5*dxy, x.z);
    vec2 index = (x-y).xy;
    
    float ch = rand(index);
    float hit = rand(index+vec2(iTime-mod(iTime,T)));
    float zap_chance = intensity;
    
    if(hit < zap_chance && ch < .6 && mod(iTime, T) < .5*T)
        sdf = add(sdf, zap(y));
    
    if(ch < .2) // capacities4
    {
        if(hit <zap_chance)
            sdf = add(sdf, vec2(cylinder((y+1.*c.yyx).yzx, vec2(.2,.6)), 7.));
        else
            sdf = add(sdf, vec2(cylinder((y+1.*c.yyx).yzx, vec2(.2,.6)), 8.));
    }
    else if(ch < .6) // chips
    {
        if(hit < zap_chance)
            sdf = add(sdf, chip(y, 7.));
        else
            sdf = add(sdf, chip(y, 6.));
    }
        
    float guard = -box(y,dxy*c.xxx);
    guard = abs(guard)+dxy*.1;
    sdf.x = min(sdf.x, guard);
    
    	    
    ////bounding sphere
    //sdf.x = max(sdf.x, sphere(x-c.yxy*t,15.));
    
    return sdf;
}

const float tmax = 100.;
const int nmax = 3000;
const float epsilon = 5.e-4;
vec2 intersect(in vec3 origin, in vec3 direction, out vec3 intersection)
{
    
    float t = 0.;
    for(int i=0; i<nmax; ++i)
    {
        intersection = origin+t*direction;
        vec2 sc = scene(intersection);
        if(sc.x<epsilon*t)return vec2(t,sc.y);
        t+=sc.x;
        if(t>tmax)return vec2(tmax, 6.);
    }
    return vec2(tmax,6.);
}

const float dx = 5.e-3;
vec3 normal(vec3 x)
{
    float sc = scene(x).x;
    return normalize(vec3(scene(x+dx*c.xyy).x-sc, scene(x+dx*c.yxy).x-sc, scene(x+dx*c.yyx).x-sc));
}

vec3 synthcol(float scale, float phase)
{
    //scale = .5+.5*sin(scale);
    vec3 c2 = vec3(207.,30.,102.)/255.,
        c3 = vec3(245., 194., 87.)/255.;
    mat3 r1 = rot((1.e-0*phase)*vec3(1.1,1.3,1.5));
    return 
        1.1*mix
        (
        	abs(r1*c2),
            abs(r1*c3), 
			scale
    	);
}

float line(vec2 x, vec2 from, vec2 to, float delta)
{
    float pm = dxy;
    
    vec2 dir = to-from, e1 = normalize(dir), e2 = vec2(e1.y, -e1.x);
    mat2 r = mat2(e1, e2), trans = inverse(r);
    vec2 p0 = trans*from, p1 = trans*(to), p = trans*x;
    
    float l = length(dir);
    
    float scale = smoothstep(p0.y+pm/2.-3.*delta, p0.y+pm/2.-delta, p.y+l)
        *(1.-smoothstep(p0.y+pm/2.+delta, p0.y+pm/2.+3.*delta, p.y+l))
        *(smoothstep(-3.*delta, -delta, p.x+.5*l))
        *(1.-smoothstep(delta, 3.*delta, p.x-.5*l));                  
    return scale;
}

float circle(vec2 x, float r)
{
    return 1.-step(length(x),r);
}

vec3 colorize(vec3 origin, vec3 intersection, float material)
{
    if(material < 0.) return c.yyy;
    vec3 color = c.yyy;
    vec3 n = normal(intersection);
    /*vec3 xi = rotate(intersection, 3.e-1*c.yyx*iTime);
    xi+=c.yxy*t;*/
    vec3 light = vec3(10.,0.,0.);
    /*light = rotate(light, 3.e-1*c.yyx*iTime);
	light+=c.yxy*t;*/
    vec3 l = normalize(light-intersection);
   
    vec3 r = normalize(reflect(-l,n));
  	vec3 v = normalize(origin-intersection);
    
    if(material <= 0.)//background
    {
        
    }
    else if(material <= 2.)//skyscraper color I
    {
        color += synthcol(5.e0*intersection.z, 62.2);
        color += synthcol(5.e0*intersection.z, 65.2)*(pow(dot(r,v),2.)+2.4*dot(l,n));
    }
    else if(material <= 3.)//skyscraper color II
    {
        color += synthcol(8.e0*intersection.z, 9.6);
        color += 2.*synthcol(8.e0*intersection.z, 14.6)*(pow(dot(r,v),2.)+2.1*dot(l,n));
    }
    else if(material <= 4.)//skyscraper color III
    {
        color += synthcol(8.e0*intersection.z, 79.64);
        color += 2.*synthcol(8.e0*intersection.z, 84.64)*(pow(dot(r,v),2.)+2.1*dot(l,n));
    }
    else if(material <= 5.)//zap
    {
        color = c.xxx;
    }
    else if(material <= 6.)//chip color
    {
        /*
        float gamma = .45;
        vec3 m = pow(abs(n),2.*c.xxx);
        vec3 x = c.xxx*mfperlin2d(xi.yz, 1., 1., 1.e2, gamma),
            y = c.xxx*mfperlin2d(xi.zx, 1., 1., 1.e2, gamma),
            z = c.xxx*mfperlin2d(xi.xy, 1., 1., 1.e2, gamma);
        color += normalize(synthcol(
            5.e-2*(pow(2.*dot(r,v),5.)+.3*dot(l,n))
        , asinh(cos(iTime))))+.2*(x*m.x+y*m.y+z*m.z)/(m.x+m.y+m.z);  
*/

        color += .2*c.xxx;
        color += .2*c.xxx*(pow(dot(r,v),6.)+10.2*dot(l,n));
//         color += .6*synthcol(5.*sin(intersection.y*intersection.x), 7.);
//         color += synthcol(5.*cos(intersection.y*intersection.x), 0.)*(pow(dot(r,v),6.)+10.2*dot(l,n));
        color *= .5;
    }
    else if(material <= 7.)//capacities
    {
        color += .2*c.xyy;
        color += .2*c.xyy*(pow(dot(r,v),6.)+10.2*dot(l,n));
//         color += .6*synthcol(5.*sin(1.e-2*intersection.y*intersection.x), iTime+3.);
//         color += .2*synthcol(5*cos(2.e-2*intersection.y*intersection.x), iTime+9)*(pow(dot(r,v),6.)+10.2*dot(l,n));
    }
    else if(material <= 8.)//house color II
    {
        color += .3*c.yxy+.1*c.xyx;
        color += .5*c.yxy*(pow(dot(r,v),6.)+10.2*dot(l,n));
    }
    else if(material <= 9.)//house color III
    {
        color += synthcol(12.e0*intersection.z, 60.3);
        color += 1.*synthcol(12.e0*intersection.z, 65.3)*(pow(dot(r,v),2.)+1.9*dot(l,n));
    }
    else if(material <= 10.)//house color IV
    {
        color += synthcol(1.e0*intersection.z, 140.71);
        color += 1.*synthcol(1.e0*intersection.z, 131.71)*(pow(dot(r,v),2.)+2.9*dot(l,n));
    }
    else if(material <= 11.)//house color V
    {
        color += synthcol(8.e0*intersection.z, 177.76);
        color += 1.*synthcol(8.e0*intersection.z, 182.76)*(pow(dot(r,v),2.)+1.9*dot(l,n));
    }
    else if(material <= 12.)//house color VI
    {
        color += synthcol(23.e0*intersection.z, 179.5);
        color += 1.*synthcol(23.e0*intersection.z, 184.5)*(pow(dot(r,v),2.)+2.9*dot(l,n));
    }
    else if(material <= 13.)//light dot
    {
        color += c.xxx*pow(dot(r,v),5.);
    }
    else if(material <= 14.) //osci lines
    {
        float pm = dxy, delta = .01;
        vec2 p = mod(intersection.xy, pm);
        float scale = max(smoothstep(pm/2.-3.*delta, pm/2.-delta, p.x)*(1.-smoothstep(pm/2.+delta, pm/2.+3.*delta, p.x)),
                          smoothstep(pm/2.-3.*delta, pm/2.-delta, p.y)*(1.-smoothstep(pm/2.+delta, pm/2.+3.*delta, p.y)));
        color += mix(c.yyy, synthcol(.5+.5*sin(3.e-1*length(intersection.xy-10.*c.yx)- 5.e0*t), iTime), scale);
    }
    else if(material <= 15.) //conductors
    {
        vec3 x = intersection;
        x += 1.*c.yxy+iTime*c.yxy;
    
    	vec3 y = vec3(mod(x.xy, dxy)-.5*dxy, x.z);
    	vec2 index = (x-y).xy;
    
        float scale = 0., delta = .04;        
        
    	float ch = rand(index);
        
    	if(ch >= .2 && ch < .6) // chips
        {
            ch = rand(index + dxy*c.xy); //right neighbor?
            if(ch >= .2 && ch < .6 || ch > .7) // chips
            {
                vec2 z = vec2(y.x, mod(y.y, delta)-.5*delta), zi=1./delta*(y.xy-z);
                if(abs(zi.y-y.y)<3.)
	                scale = max(scale, line(z-.5*dxy, dxy*c.xy, c.yy, .005));
            }
            
            ch = rand(index - dxy*c.xy); //left neighbor?
            if(ch >= .2 && ch < .6 || ch > .7) // chips
            {
                vec2 z = vec2(y.x, mod(y.y, delta)-.5*delta), zi=1./delta*(y.xy-z);
                if(abs(zi.y-y.y)<3.)
	                scale = max(scale, line(z+.5*dxy, -dxy*c.xy, c.yy, .005));
            }
            
            ch = rand(index + dxy*c.yx); //upper neighbor?
            if(ch >= .2 && ch < .6 || ch > .7) // chips
            {
                vec2 z = vec2(mod(y.x, delta)-.5*delta, y.y), zi=1./delta*(y.xy-z);
                if(abs(zi.x-y.x)<3.)
	                scale = max(scale, line(z-.5*dxy, -dxy*c.yx, c.yy, .005));
            }
            
            ch = rand(index - dxy*c.yx); //lower neighbor?
            if(ch >= .2 && ch < .6 || ch > .7) // chips
            {
                vec2 z = vec2(mod(y.x, delta)-.5*delta, y.y), zi=1./delta*(y.xy-z);
                if(abs(zi.x-y.x)<3.)
	                scale = max(scale, line(z+.5*dxy, dxy*c.yx, c.yy, .005));
            }
        }
        
        
        
        color = scale*c.yxy*(.75+.25*sin(5.*iTime));
        
        /*
        mat2 ra = mat2(cos(iTime), sin(iTime), -sin(iTime), cos(iTime));
        float scale = max(line(intersection.xy, ra*c.xy, ra*c.xx),
                          line(intersection.xy, ra*c.xx, ra*c.yx));
        scale = max(scale, line(intersection.xy, ra*c.xy, ra*c.yy));
        scale = max(scale, line(intersection.xy, ra*c.yy, ra*c.yx));
        color = scale*c.xyy;*/
    }
    float d = length(origin-intersection);
    color = mix(color, mix(c.xyy,  c.yyx, .5+.3*sin(6.*iTime)), clamp(1.-d,0.,1.));
    
    //fog
    /*
    color = mix(color, mix(c.yxy,  c.yyx, .5+.3*sin(1.*iTime-.5*intersection.y)), clamp(.02*d*d,0.,1.));
    color += mix(color, mix(c.xyy,  c.yxy, .5+.3*sin(2.*iTime-.5*intersection.y)), clamp(.02*d*d,0.,1.));
    color += mix(color, mix(c.yyx,  c.xyy, .5+.3*sin(3.*iTime-.5*intersection.y)), clamp(.02*d*d,0.,1.));
    color *= .3;
    */
    return color;
}

vec3 render(float time, vec2 fragCoord)
{
    t = time;
 
    //vec3 ray_origin = c.yyx*(.25-scene(.4*c.yyx).x);
    /*vec3 ray_origin = c.yyx*(.25-scene(.4*c.yyx).x);
    vec3 x = rotate(ray_origin, 3.e-1*c.yyx*iTime);
	x+=c.yxy*t;
    float p0 = mfperlin2d(x.xy, 3., 1., 21.05, .15);
    ray_origin = c.yyx*(.25+p0);*/
    vec3 ray_origin = c.yyx-c.yxy;
    vec3 ray_target = c.yyy;
    vec3 camera_right = c.xyy;
    vec3 camera_up = normalize(cross(ray_target-ray_origin, -camera_right));
    
    vec3 color=c.yyy;
    
#if AA!=1
    for(int i=0; i<AA; ++i)
    	for(int j=0; j<AA; ++j)
        {
            vec2 o = vec2(float(i),float(j)) / float(AA) - 0.5;
        	vec2 p = (-iResolution.xy + 2.*(fragCoord+o))/iResolution.y;
#else 
            vec2 p = (-iResolution.xy + 2.*fragCoord)/iResolution.y;
    
#endif
    		vec3 ray_direction = normalize(ray_target+p.x*camera_right+p.y*camera_up-ray_origin);        
            
            vec3 intersection;
            vec2 mat_t = intersect(ray_origin, ray_direction, intersection);
            
            color +=colorize(ray_origin, intersection, mat_t.y);
#if AA!=1
        }
    color/=float(AA*AA);
#else
#endif
    
    return color;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    
    //motion blur
    //vec3 color = .5*(render(iTime, fragCoord)+render(iTime+.1+.1*sin(2.*iTime), fragCoord));

    //simple
    vec3 color = render(iTime, fragCoord);
    
    fragColor = vec4(color, 1.);
}



void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
