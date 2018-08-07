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

#define AA 1

#define T .5

const vec2 c = vec2(1.,0.);
const float pi = 3.14159;
float t;

uniform float iTime;
uniform vec2 iResolution;
uniform float intensity;

float rand(vec2 a0)
{
    return fract(sin(dot(a0.xy ,vec2(12.9898,78.233)))*43758.5453);
}

vec3 taylorInvSqrt(vec3 r) 
{     
    return 1.79284291400159-0.85373472095314*r; 
}
vec3 permute(vec3 x)
{
    return mod((x*34.+1.)*x, 289.);
}

float snoise(vec2 P) 
{     
    const vec2 C = vec2 (0.211324865405187134, 0.366025403784438597);  
    vec2 i = floor(P+dot(P, C.yy)) ; 
    vec2 x0 = P-i+dot(i, C.xx) ; 
    // Other  corners 
    vec2 i1 ; 
    i1.x = step ( x0.y , x0.x ) ;  //  1.0  i f  x0 . x > x0 . y ,  e l s e  0.0 
    i1.y = 1.0 - i1.x ; 
    // x1 = x0 − i1 + 1.0 ∗ C. xx ;  x2 = x0 − 1.0 + 2.0 ∗ C. xx ; 
    vec4 x12 = x0.xyxy + vec4 ( C.xx , C.xx * 2.0 - 1.0) ; 
    x12.xy -= i1 ; 
    //  Permutations 
    i = mod( i ,  289.0) ;  // Avoid  truncation  in  polynomial  evaluation 
    vec3 p = permute ( permute ( i.y + vec3 (0.0 , i1.y ,  1.0  ) ) + i.x + vec3 (0.0 , i1.x ,  1.0  ) ) ; 
    //  Circularly  symmetric  blending  kernel
    vec3 m = max(0.5 - vec3 ( dot ( x0 , x0 ) ,  dot ( x12.xy , x12.xy ) , dot ( x12.zw , x12.zw ) ) ,  0.0) ; 
    m = m * m ; 
    m = m * m ; 
    //  Gradients  from 41  points  on a  line ,  mapped onto a diamond 
    vec3 x = fract ( p * (1.0  /  41.0) ) * 2.0 - 1.0  ; 
    vec3 gy = abs ( x ) - 0.5  ; 
    vec3 ox = floor ( x + 0.5) ;  // round (x)  i s  a GLSL 1.30  feature 
    vec3 gx = x - ox ; //  Normalise  gradients  i m p l i c i t l y  by  s c a l i n g m 
    m *= taylorInvSqrt ( gx * gx + gy * gy ) ; // Compute  f i n a l  noise  value  at P 
    vec3 g ; 
    g.x = gx.x * x0.x + gy.x * x0.y ; 
    g.yz = gx.yz * x12.xz + gy.yz * x12.yw ; 
    //  Scale  output  to  span  range  [ − 1 ,1] 
    //  ( s c a l i n g  f a c t o r  determined by  experiments ) 
    return  130.0 * dot ( m , g ) ; 
}

float mfsimplex2d(vec2 x, float fmin, float fmax, float phi)
{
    float sum = 0.;
    float a = 1.;
    
    for(float f = fmin; f<fmax; f = f*2.)
    {
        sum = a*snoise(f*x) + sum;
        a = a*phi;
    }
    
    return sum;
}

float g(float a, float a0)
{
    return (2.*rand(vec2(a,a0))-1.);
}

mat3 Rx(float x)
{
    return mat3(c.xyyy, cos(x), sin(x), 0., -sin(x), cos(x));
}

mat3 Ry(float x)
{
    return mat3(cos(x), 0., -sin(x), c.yxy, sin(x), 0., cos(x));
}

mat3 Rz(float x)
{
    return mat3(cos(x), sin(x), 0., -sin(x), cos(x), c.yyyx);
}

vec3 rotate(vec3 a, vec3 a0)
{
    return Rz(-(a0).z)*Ry(-(a0).y)*Rx(-(a0).x)*(a);
}

float sphere(vec3 a, float a0)
{
    return length(a)-a0;
}

float box(vec3 a, vec3 a0)
{
    return length(max(abs(a)-(a0),0.));
}

vec2 box2(vec3 a, vec3 a0, float a1) 
{
    return vec2(box(a,a0),a1);
}
vec2 plane2(vec3 a, float a0, float a1)
{
    return vec2((a).z-(a0),a1);
}

vec2 add(vec2 a, vec2 a0) 
{
    return mix(a0,a,step((a).x,(a0).x));
}
vec2 sub(vec2 a, vec2 a0)
{
    return mix(a0,a,step((a).x,-(a0).x));
}


float skyscraper1(vec3 x, vec3 a0, float height)
{
    
    vec3 y = vec3(x.xy, mod(x.z,a0.z)-.5*a0.z);
    
    //floors
    float sdf = box(y,vec3(a0.xy, .25*a0.z));
    
    //guard for floors
    
    float guard = -box(y,vec3(a0.xy, a0.z));
    guard = abs(guard)+.001;
    sdf = min(sdf, guard);
    
    //center
    sdf = min(sdf, box(x,vec3(.8*a0.xy, height)));
    sdf = max(sdf, box(x,vec3(.9*a0.xy, height)));
    
    return sdf;
}

vec2 scene(vec3 x)
{
    vec2 sdf = c.xy;
    
    //sdf = add(sdf, vec2(sphere(x-.5*vec3(0.,2.5,1.), .1), 13.));
    
    x = rotate(x, 3.e-1*c.yyx*iTime+.01*sin(iTime)*c.xyy);
	x+=c.yxy*t;
	
    //ground
    //float p0 = mfperlin2d(x.xy, 3., 1., 21.05, .45);
    float p0 = mfsimplex2d(x.xy, .24, 82.05, .35+.1*sin(iTime/T));
    sdf = add(sdf, plane2(x, p0, 6.));
   
    //skyscrapers
    vec3 z = vec3(mod(x.xy, .1)-.05, x.z);
    vec2 index = x.xy-z.xy;
    //float p1 = mfperlin2d(index.xy, 3., 1., 21.05, .45);
    float p1 = mfsimplex2d(index.xy, .24, 82.05, .35+.1*sin(iTime/T));
    float rotation = rand(index)*2.*pi;
    vec3 y = rotate(z,rotation*c.yyx);
    float dh = .04*sin(4.*t/T-10.*length(index))+.2*rand(index);
    //dh = mod(dh, .4);
    float height = (.55+.35*sin(p1)+dh+(.5+.5*rand(index))*intensity);
    float width = 2.*( .014+.004*rand(2.*index));
    float color = 2.+floor(3.*rand(index));
    if(p1 < -.12)
	    sdf = add(sdf, vec2(skyscraper1(y-p1*c.yyx, vec3(width*c.xx,height/52.), height), color));
    if(p1 > -.12 && p1 < -.02)
	    sdf = add(sdf, vec2(skyscraper1(y-p1*c.yyx, vec3(width*c.xx,height/52.), .8*height), color+5.));
    if(p1 > 0. && p1 < .12)
    //else
        sdf = add(sdf, vec2(skyscraper1(y-p1*c.yyx, vec3(width*c.xx,height/52.), .4*height), color+8.));
    
	float guard = -box2(z-p1*c.yyx, vec3(.1,.1,10.*height), -1.).x;
    guard = abs(guard)+1.*0.01;
    sdf.x = min(sdf.x, guard);
        
//*/
    
    
    //bounding sphere
    sdf.x = max(sdf.x, sphere(x-c.yxy*t,10.));
    
    return sdf;
}

const float tmax = 10.;
const int nmax = 3000;
const float epsilon = 1.e-4;
vec2 intersect(in vec3 origin, in vec3 direction, out vec3 intersection)
{
    
    float t = 0.;
    for(int i=0; i<nmax; ++i)
    {
        if(t<0.)return vec2(tmax,-1.);
        intersection = origin+t*direction;
        vec2 sc = scene(intersection);
        if(sc.x<epsilon*t)return vec2(t,sc.y);
        t+=sc.x;
        if(t>tmax)return vec2(tmax, -1.);
    }
    return vec2(tmax,-1.);
}

const float dx = 2.5e-5;
vec3 normal(vec3 x)
{
    float sc = scene(x).x;
    return normalize(vec3(scene(x+dx*c.xyy).x-sc, scene(x+dx*c.yxy).x-sc, scene(x+dx*c.yyx).x-sc));
}

mat3 rot(vec3 p)
{
    vec3 cp = cos(p), sp = sin(p);
    mat3 m = mat3(cp.y*cp.x, cp.x*sp.z+cp.z*sp.x*sp.y, sp.x*sp.z-cp.x*cp.z*sp.y, 
           -cp.y*sp.z, cp.x*cp.z-sp.x*sp.y*sp.z, cp.z*sp.x+cp.x*sp.y*sp.z, 
           sp.y, -cp.y*sp.x, cp.x*cp.y);
    return m;
}

vec3 synthcol(float scale, float phase)
{
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

const float alpha = 8.;
vec3 colorize(vec3 origin, vec3 intersection, float material)
{
    //if(material < 0.) return .2*c.xxx;
    vec3 color = c.yyy;
    vec3 n = normal(intersection);
    vec3 light = 20.7*c.yyx;
    /*
    vec3 xi = rotate(intersection, 3.e-1*c.yxx*iTime);
    xi+=c.yxy*t;
    
    light = rotate(light, 3.e-1*c.yyx*iTime);
	light+=c.yxy*t;*/
    vec3 l = normalize(light-intersection);
   
    //vec3 r = normalize(reflect(-1.e2*l,1.e2*n));
  	vec3 v = normalize(origin-intersection);
    
    vec3 r = normalize(l+v);
    if(material <= 0.)//background
    {
        color += .5*c.xxx;
    }
    else if(material <= 2.)//skyscraper color I
    {
        color += synthcol(5.e0*intersection.z, 62.2);
        color += synthcol(5.e0*intersection.z, 65.2)*(pow(dot(r,v),alpha)+2.4*dot(l,n));
    }
    else if(material <= 3.)//skyscraper color II
    {
        color += synthcol(5.e0*intersection.z, 9.6);
        color += 2.*synthcol(5.e0*intersection.z, 14.6)*(pow(dot(r,v),alpha)+2.1*dot(l,n));
    }
    else if(material <= 4.)//skyscraper color III
    {
        color += synthcol(5.e0*intersection.z, 79.64);
        color += 2.*synthcol(5.e0*intersection.z, 84.64)*(pow(dot(r,v),alpha)+2.1*dot(l,n));
    }
    else if(material <= 6.)//mountain color
    {
        
        /*float gamma = .45;
        vec3 m = pow(abs(n),2.*c.xxx);
        vec3 x = c.xxx*mfperlin2d((origin-intersection).yz, 1., 1., 1.e2, gamma),
            y = c.xxx*mfperlin2d((origin-intersection).zx, 1., 1., 1.e2, gamma),
            z = c.xxx*mfperlin2d((origin-intersection).xy, 1., 1., 1.e2, gamma);
        color -= 1.*normalize(synthcol(
            5.e-2*(pow(2.*dot(r,v),5.)+.3*dot(l,n))
        , asinh(cos(iTime))))+.2*(x*m.x+y*m.y+z*m.z)/(m.x+m.y+m.z);  */

        color += 1.e-1*synthcol(1.e-2*length(origin-intersection), .5+.5*sin(2.*iTime));
        color += synthcol(1.e-2*length(origin-intersection), .5+.5*sin(iTime))*(5.e0*pow(dot(n,v), alpha)+1.3e-2*dot(l,n));
        color *= .5;
        //color += 1.e-1*synthcol(1.e0*(1.e-2*(origin-intersection).z), iTime)*dot(l,n);
        //color += 1.e-1*(synthcol(1.e0*(1.e-2*(origin-intersection).z), iTime+5.))*pow(dot(r,v),alpha)*65.;
         // 10.2));
    }
    else if(material <= 7.)//house color I
    {
        color += synthcol(15.e0*intersection.z, 29.14);
        color += 1.*synthcol(5.e0*intersection.z, 34.14)*(12.*pow(dot(r,v),alpha)+11.9*dot(l,n));
    }
    else if(material <= 8.)//house color II
    {
        color += synthcol(15.e0*intersection.z, 94.22);
        color += 1.*synthcol(5.e0*intersection.z, 99.22)*(pow(dot(r,v),alpha)+1.9*dot(l,n));
    }
    else if(material <= 9.)//house color III
    {
        color += synthcol(15.e0*intersection.z, 60.3);
        color += 1.*synthcol(5.e0*intersection.z, 65.3)*(pow(dot(r,v),alpha)+1.9*dot(l,n));
    }
    else if(material <= 10.)//house color IV
    {
        color += synthcol(15.e0*intersection.z, 140.71);
        color += 1.*synthcol(5.e0*intersection.z, 131.71)*(pow(dot(r,v),alpha)+2.9*dot(l,n));
    }
    else if(material <= 11.)//house color V
    {
        color += synthcol(15.e0*intersection.z, 177.76);
        color += 1.*synthcol(5.e0*intersection.z, 182.76)*(pow(dot(r,v),alpha)+1.9*dot(l,n));
    }
    else if(material <= 12.)//house color VI
    {
        color += synthcol(15.e0*intersection.z, 179.5);
        color += 1.*synthcol(5.e0*intersection.z, 184.5)*(pow(dot(r,v),alpha)+2.9*dot(l,n));
    }
    else if(material <= 13.)//light dot
    {
        color += c.xxx*pow(dot(r,v),5.);
    }
    float d = length(origin-intersection);
    color = mix(color, mix(c.xyy,  c.yyx, .5+.3*sin(6.*iTime)), clamp(1.-d,0.,1.));
    
    //fog
    
    d += 1.;
    color = mix(color, mix(synthcol(.04*d*d,iTime+30.),  synthcol(.04*d*d,iTime+40.), .5+.3*sin(1.*iTime-1.5e-1*length(intersection))), clamp(.04*d*d,0.,1.));
    color += mix(color, mix(synthcol(.04*d*d,iTime+20.),  synthcol(.04*d*d,iTime+50.), .5+.3*sin(2.*iTime-2.5e-1*length(intersection))), clamp(.04*d*d,0.,1.));
    color += mix(color, mix(synthcol(.04*d*d,iTime+10.),  synthcol(.04*d*d,iTime+60.), .5+.3*sin(3.*iTime-3.5e-1*length(intersection))), clamp(.04*d*d,0.,1.));
    color *= .2;
    color = clamp(color, 0., 1.);
    return color;
}

vec3 render(float time, vec2 fragCoord)
{
    //+90
    t = .5*time+1.-40.;
 
    //vec3 ray_origin = c.yyx*(.25-scene(.4*c.yyx).x);
    vec3 ray_origin = c.yyx*(.65-scene(.4*c.yyx).x);
    vec3 x = rotate(ray_origin, 3.e-1*c.yyx*iTime);
	x+=c.yxy*t;
    //float p0 = mfperlin2d(x.xy, 3., 1., 21.05, .15);
    float p0 = mfsimplex2d(x.xy, .24, 42.05, .15);
    ray_origin = c.yyx*(.65+p0);
    vec3 ray_target = ray_origin + 1.4*c.yxy-1.*c.yyx;
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
