const char *gfx_source = "#version 130\nuniform float a;vec2 b,c,d=vec2 (1280.,720.);vec3 e=vec3 (1.,0.,-1.),f,g,h,i,j,k;float l=0.,m=acos(-1.);int n=550,o,p;float q(vec2 r){vec3 s=fract(vec3 (r.xyx)*.1031);s+=dot(s,s.yzx+33.33);return fract((s.x+s.y)*s.z);}float t(vec2 u){vec2 o=floor(u);u=fract(u);u=smoothstep(e.yy,e.xx,u);vec2 v=vec2 (q(o),q(o+e.xy)),w=vec2 (q(o+e.yx),q(o+e.xx));v=e.zz+2.*mix(v,w,u.y);return mix(v.x,v.y,u.x);}float x(vec2 h,float l,float y,float z){float i=0.,A=1.,B=0.,C;for(float D=l;D<y;D*=2.){C=t(D*h);i+=A*C;A*=z;B+=1.;}i*=(1.-z)/(1.-pow(z,B));return i;}vec2 E(vec2 F,vec2 G){return (F.x<G.x)?F:G;}vec2 H(in vec3 I){I.y+=.7*a;float J=x(I.xy,4.,450.,mix(.25,.5,smoothstep(.05,.35,I.z))),K=J,L;J=.5+.5*J;J*=smoothstep(.7,-.1,abs(I.x));J+=.1*K;vec2 M=vec2 (I.z-.6*J,1.),N=mod(I.xy,.4)-.2,O=I.xy-N;L=q(O);if(L>.9){M=E(M,vec2 (abs(length(vec3 (N,I.z-L+.9))-2.*(L-.9))-.001,2.));}vec2 A=vec2 (mod(I.x,.05)-.025,I.z);M=E(M,vec2 (length(A)-.005,3.));M=E(M,vec2 (abs(length(A)-.025)-.001,4.));return M;}vec3 P(in vec3 I){float Q=1.e-4,c=H(I).x;return normalize(vec3 (H(I+Q*e.xyy).x-c,H(I+Q*e.yxy).x-c,H(I+Q*e.yyx).x-c));}vec3 R(){float S=dot(k,i),T=abs(dot(reflect(normalize(h-k),i),g));if(c.y==1.){return mix(.7*vec3 (0.38,0.39,0.61),vec3 (0.22,0.01,0.20),smoothstep(.0,.5,h.z))+.1*vec3 (0.56,0.12,0.47)*S+e.xxx*pow(T,6.);}else if(c.y==2.){return .2*e.yyx+.1*e.xxx*S+8.*vec3 (0.78,0.39,0.31)*pow(T,2.);}else if(c.y==3.){return .1*vec3 (1.00,0.87,0.74)+.1*vec3 (1.00,0.87,0.74)*S+2.*vec3 (1.00,0.87,0.74)*pow(T,2.);}else if(c.y==4.){return e.xyy+e.xyy*S+e.xyy*pow(T,2.);}}void main(){b=(gl_FragCoord.xy-.5*d.xy)/d.yy;f=vec3 (b,1.);g=e.yyz;float r=m/4.;vec2 U=vec2 (cos(r),sin(r));mat2 V=mat2 (U.x,U.y,-U.y,U.x);g.yz=V*g.yz;for(o=0;o<n;++o){h=f+l*g;c=H(h);if(c.x<1.e-4)break;l+=min(c.x,1.e-2);}if(o<n){i=P(h);k=h+.1*i;j=R();if(c.y==2.||c.y==4.){for(p=0;p<2;++p){f=h;g=refract(g,i,.5);l=.1;for(o=0;o<n;++o){h=f+l*g;c=H(h);if(c.x<1.e-4)break;l+=min(c.x,5.e-2);}if(o<n){i=P(h);k=h+.1*i;vec3 W=R();j=mix(j,W,.7);}}}}j*=j;gl_FragColor=vec4 (clamp(j,0.,1.),1.0);}\0";