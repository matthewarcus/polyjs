#if 0
#version 300 es
precision highp float;
uniform float iTime;
uniform vec2 iResolution;
uniform vec4 iMouse;

out vec4 outColor;

void mainImage( out vec4 fragColor, in vec2 fragCoord );

void main() {
  mainImage(outColor, gl_FragCoord.xy);
}
#endif

////////////////////////////////////////////////////////////////////////////////
//
// Goursat Quartic Surfaces
//
// http://mathworld.wolfram.com/GoursatsSurface.html
// https://www.mathcurve.com/surfaces/goursat/goursat.shtml
//
// Controls:
// <mouse>: change view direction
// <up/down>: zoom
// <left/right>: cycle through some interesting parameters
// 'c': clip at z = 0
// 'g': show grid
// 'r': do rotation
//
// Quartic surfaces with octahedral symmetry. Surface (including normals)
// is raytraced using analytic solution to quartic due to Lanczos and Kahan.
//
////////////////////////////////////////////////////////////////////////////////

uniform mat4 iView; // Not used for Shadertoy - view matrix in VR

// PARAMS
bool dorotate = true;

// Lighting
vec3 light = normalize(vec3(1,1,1));
float ambient = 0.4;
float diffuse = 0.6;
float specular = 0.8;
float specularpow = 10.0;
vec3 specularcolor = vec3(1);

// Debug
bool alert = false;

void assert(bool t) {
  if (!t) alert = true;
}

bool eq(float x, float y) {
  return abs(x-y) < 1e-4;
}

bool eq(vec4 p, vec4 q) {
  return eq(p.x,q.x) && eq(p.y,q.y) && eq(p.z,q.z) && eq(p.w,q.w);
}

bool eq(mat4 m, mat4 n) {
  return eq(m[0],n[0]) && eq(m[1],n[1]) && eq(m[2],n[2]) && eq(m[3],n[3]);
}

const float PI =  3.141592654;

vec2 rotate(vec2 p, float t) {
  return p * cos(t) + vec2(p.y, -p.x) * sin(t);
}

// The Kahan algorithm, for explanation see:
// https://people.eecs.berkeley.edu/~wkahan/Math128/Cubic.pdf

float sgn(float x) {
  return x < 0.0? -1.0: 1.0; // Return 1 for x == 0
}

int quadratic(float A, float B, float C, out vec2 res) {
  float x1,x2;
  float b = -0.5*B;
  float q = b*b - A*C;
  if (q < 0.0) return 0;
  float r = b + sgn(b)*sqrt(q);
  if (r == 0.0) {
    x1 = C/A; x2 = -x1;
  } else {
    x1 = C/r; x2 = r/A;
  }
  res = vec2(x1,x2);
  return 2;
}

int quadratic(vec3 coeffs, out vec2 res) {
  return quadratic(coeffs[0],coeffs[1],coeffs[2],res);
}
  
void eval(float X, float A, float B, float C, float D,
          out float Q, out float Q1, out float B1,out float C2) {
  float q0 = A*X;
  B1 = q0+B;
  C2 = B1*X+C;
  Q1 = (q0+B1)*X + C2;
  Q = C2*X + D;
}

// Solve: Ax^3 + Bx^2 + Cx + D == 0
// Find one real root, then reduce to quadratic.
int cubic(float A, float B, float C, float D, out vec3 res) {
  float X,b1,c2;
  if (A == 0.0) {
    X = 1e8; A = B; b1 = C; c2 = D;
  } else if (D == 0.0) {
    X = 0.0; b1 = B; c2 = C;
  } else {
    X = -(B/A)/3.0;
    float t,r,s,q,dq,x0;
    eval(X,A,B,C,D,q,dq,b1,c2);
    t = q/A; r = pow(abs(t),1.0/3.0); s = sgn(t);
    t = -dq/A; if (t > 0.0) r = 1.324718*max(r,sqrt(t));
    x0 = X - s*r;
    if (x0 != X) {
      // Newton-Raphson
      for (int i = 0; i < 6; i++) {
        X = x0;
        eval(X,A,B,C,D,q,dq,b1,c2);
        if (dq == 0.0) break;
        x0 -= (q/dq);
      }
      if (abs(A)*X*X > abs(D/X)) {
        c2 = -D/X; b1 = (c2 - C)/X;
      }
    }
  }
  res.x = X;
  return 1 + quadratic(A,b1,c2,res.yz);
}

int cubic(vec4 coeffs, out vec3 res) {
  float A = coeffs[0], B = coeffs[1], C = coeffs[2], D = coeffs[3];
  return cubic(A,B,C,D,res);
}

// Special wrapper for cubic function for solving quartic.
// Find largest real root of x**3 + a*x**2 + b*x + c
// Assume c < 0
float qcubic(float a, float b, float c) {
  // c is always <= 0, but may be very
  // small, in which case we return an
  // approximation. Never return < 0.
  assert(c <= 0.0);
  if (c == 0.0) return 0.0;
  if (false) {
    // This helps with double roots, but sometimes is
    // completely the wrong thing to do.
    // Further investigation required.
    if (c > -1e-6) {
      if (b > 1e-10) return -c/b;
      //if (b > 0.0) return -c/b; // Keep it simple.
      if (b > -1e-4) return 0.0;
    }
  }
  vec3 res;
  int nroots = cubic(1.0,a,b,c,res);
  if (nroots == 1) return res.x;
  else return max(res.x,max(res.y,res.z));
}

int quartic(vec4 coeffs, out vec4 res) {
  float c1 = coeffs[0];
  float c2 = coeffs[1];
  float c3 = coeffs[2];
  float c4 = coeffs[3];
  float alpha = 0.5*c1;
  float A = c2-alpha*alpha;
  float B = c3-alpha*A;
  float a,b,beta,psi;
  psi = qcubic(2.0*A-alpha*alpha, A*A+2.0*B*alpha-4.0*c4, -B*B);
  //assert(!isnan(psi));
  //assert(!isinf(psi));
  //assert(psi >= 0.0);
  a = sqrt(max(0.0,psi));
  beta = 0.5*(A + psi);
  if (psi <= 0.0) {
    b = sqrt(max(beta*beta-c4,0.0));
  } else {
    b = 0.5*a*(alpha-B/psi);
  }
  int resn = quadratic(1.0,alpha+a,beta+b,res.xy);
  vec2 tmp;
  if (quadratic(1.0,alpha-a,beta-b,tmp) != 0) { 
    res.zw = res.xy;
    res.xy = tmp;
    resn += 2;
  }
  return resn;
}

int quartic(float A, float B, float C, float D, float E, out vec4 roots) {
  int nroots;
  // There may be a better heuristic for this.
  // but this avoids the worst glitches.
  if (abs(E) < 10.0*abs(A)) {
    vec4 coeffs = vec4(B,C,D,E)/A;
    nroots = quartic(coeffs,roots);
  } else {
    // It can be advantageous to use the coefficients in the
    // opposite order, thus solving for the reciprocal.
    vec4 coeffs = vec4(D,C,B,A)/E;
    nroots = quartic(coeffs,roots);
    for (int i = 0; i < 4; i++) {
      if (i == nroots) break;
      roots[i] = 1.0/roots[i];
    }
  }
  return nroots;
}

struct Surface {
  vec4 params;
  vec3 p;
  int colorscheme;
};

// Equation: pp.pp + k(p.p)^2 + k'a^2(p.p) + k''a^4 = 0
// Derivative: 4ppp + 4k(p.p)p + 2k'a^2p
// Expansion with p => p+tr:
// pp => (p+tr)(p+tr) = pp + 2tpr + t^2rr
// pp.pp => (pp + 2tpr + t^2rr).(pp + 2tpr + t^2rr)
//  = pp.pp + 4tpp.pr + 6t^2pp.rr + 4t^3pr.rr + t^4rr.rr 
// p.p  => (p+tr).(p+tr) = p.p + 2tp.r + t^2r.r = p.p + 2tp.r + t^2
// (p.p)^2 = (p.p + 2tp.r + t^2)(p.p + 2tp.r + t^2) =
//         = p.p^2 + 4t^2(p.r)^2 + t^4 + 2(2t(p.p)(p.r) + (p.p)t^2 + 2t^3(p.r))
//         = p.p^2 + 4t^2(p.r)^2 + t^4 + 4t(p.p)(p.r) + 2(p.p)t^2 + 4t^3(p.r))
// ie.
// pp.pp + 4tpp.pr + 6t^2pp.rr + 4t^3pr.rr + t^4rr.rr +
// k(p.p^2 + t4[(p.p)(p.r)] + t^2[4(p.r)^2 + 2(p.p)] + t^3[4(p.r)] + t^4) +
// k'a^2(p.p + 2tp.r + t^2) +
// k''a^4
// collecting terms:
// t^0: pp.pp +   k(p.p)^2 +             k'a^2(p.p) + k''a^4
// t^1: 4pp.pr + 4k(p.p)(p.r) +         2k'a^2(p.r)
// t^2: 6pp.rr +  k[4(p.r)^2 + 2(p.p)] + k'a^2
// t^3: 4pr.rr + 4k(p.r)
// t^4: rr.rr +   k
//
// Can adjust p so that p.r = 0 & that simplifies things.

int goursatsurface(Surface surface, vec3 p, vec3 r, out vec4 roots) {
  float k = surface.params[0];
  float k1 = surface.params[1];
  float k2 = surface.params[2];
  float a = surface.params[3];
  vec3 pp = p*p;
  vec3 pr = p*r;
  vec3 rr = r*r;
  float p2 = dot(p,p);
  float a2 = a*a;
  float a4 = a2*a2;
  // Check we have adjusted p so that p.r = 0!
  //assert(eq(dot(p,r),0.0));
  //assert(eq(dot(r,r),1.0));
#if 0
  // Unoptimized version
  float pdr = dot(p,r);
  float A = dot(rr,rr) + k;
  float B = 4.0*dot(pr,rr) + 4.0*k*pdr;
  float C = 6.0*dot(pp,rr) +     k*(4.0*pdr*pdr + 2.0*p2) + k1*a2;
  float D = 4.0*dot(pp,pr) + 4.0*k*(p2*pdr)           + 2.0*k1*a2*pdr;
  float E = dot(pp,pp)         + k*p2*p2 + k1*a2*p2 + k2*a4;
#else
  // Optimized, with p.r = 0
  float A =     dot(rr,rr) +     k;
  float B = 4.0*dot(pr,rr);
  float C = 6.0*dot(pp,rr) + 2.0*k*p2    + k1*a2;
  float D = 4.0*dot(pp,pr);
  float E =     dot(pp,pp) +     k*p2*p2 + k1*a2*p2 + k2*a4;
#endif
  return quartic(A,B,C,D,E,roots);
}

vec3 goursatnormal(Surface surface, vec3 p) {
  // 4ppp + 4k(p.p)p + 2k'a^2p
  float k = surface.params[0];
  float k1 = surface.params[1];
  float a = surface.params[3];
  return 4.0*p*p*p + 4.0*k*dot(p,p)*p + 2.0*k1*a*a*p;
}

vec3 applylighting(vec3 baseColor, vec3 p, vec3 n, vec3 r) {
  if (dot(r,n) > 0.0) n = -n; // Face forwards
  vec3 c = baseColor*ambient;
  c += baseColor*diffuse*(max(0.0,dot(light,n)));
  float s = pow(max(0.0,dot(reflect(light,n),r)),specularpow);
  c += specular*s*specularcolor;
  return c;
}

struct Result {
  vec3 p;
  vec3 n;
  vec3 basecolor;
  float t;
};

float gridline(vec3 p) {
  // Draw some gridlines on surface
  vec3 t = fract(p*4.0);
  t = min(t,1.0-t);
  float d = min(t.x,min(t.y,t.z));
  return smoothstep(0.02,0.025,d);
}

int dosurface(Surface surface, vec3 p0, vec3 r, out vec4 roots) {
  return goursatsurface(surface,p0,r,roots);
}
  
vec3 getnormal(Surface surface, vec3 p) {
  return goursatnormal(surface,p);
}
  
bool solve(Surface surface, vec3 p, vec3 r, float tmin, inout Result result) {
  vec4 roots;
  int nroots = dosurface(surface,p,r,roots);
  // Find smallest root greater than tmin.
  float t = result.t;
  for (int i = 0; i < 4; i++) {
    if (i == nroots) break;
    if (roots[i] > tmin && roots[i] < t) {
      t = roots[i];
    }
  }
  if (t == result.t) return false;
  p += t*r;
  vec3 n = getnormal(surface, p);
  if (dot(n,r) > 0.0) n = -n;
  n = normalize(n);
  vec3 basecolor = abs(n);
  if (surface.colorscheme == 1) {
    basecolor *= gridline(p);
  }
  result = Result(p,n,basecolor,t);
  return true;
}

// Interesting parameters from:
// https://www.mathcurve.com/surfaces.gb/goursat/goursat.shtml
vec4 goursatparams(int i) {
  int index = 0;
  if (i == index++) return vec4(0,-1,0,1);
  if (i == index++) return vec4(0,-2,2,1);
  if (i == index++) return vec4(-1,1,1,1);
  if (i == index++) return vec4(-1,-0.25,0.25,1);
  if (i == index++) return vec4(-0.5,-1,0.5,1);
  if (i == index++) return vec4(-0.5,1,-1.5,1);
  if (i == index++) return vec4(-1,4,-6,1);
  if (i == index++) return vec4(-1,1,1,1);
  if (i == index++) return vec4(-1,2,-2,1);
  else return vec4(-0.333,-0.666,0.666,1);
}

int nparams = 10;

int imod(int n, int m) {
    return n-n/m*m;
}

bool scene(vec3 p0, vec3 r, out vec3 col) {
  // Solve from closest point to origin.
  // This make p0.r = 0.
  float tmin = -dot(p0,r);
  p0 += tmin*r;
  Result res = Result(vec3(0),vec3(0),vec3(0),1e8);
  float ttime = 0.5*iTime;
  float rtime0 = floor(ttime);
  float rtime1 = floor(ttime+0.5);
  ttime = fract(2.0*ttime);
  vec4 params = vec4(0);
  params = mix(goursatparams(int(rtime0)%nparams),
               goursatparams(int(rtime1)%nparams),
               ttime);
  Surface surface = Surface(params,vec3(0),1);
  if (!solve(surface,p0,r,-tmin,res)) return false;
  col = applylighting(res.basecolor,res.p,res.n,r);
  return true;
}

vec3 transform(in vec3 p) {
#if 0
  if (iMouse.x > 0.0) {
    float theta = (2.0*iMouse.y-iResolution.y)/iResolution.y*PI;
    float phi = (2.0*iMouse.x-iResolution.x)/iResolution.x*PI;
    p.yz = rotate(p.yz,theta);
    p.zx = rotate(p.zx,-phi);
  }
#endif
  if (dorotate) {
    float t = iTime;
    p.yz = rotate(p.yz, 0.1*t);
    p.zx = rotate(p.zx, 0.222*t);
  }
  return p;
}

void mainFun(out vec4 fragColor, vec3 p, vec3 r, vec3 rcentre) {
  p = transform(p);
  r = normalize(transform(r));
  // rcentre should be ray through centre of clipspace.
  rcentre = normalize(transform(rcentre));
  // Screenspace ray derivatives
  vec3 drdx = dFdx(r);
  vec3 drdy = dFdy(r);
  vec3 color = vec3(0);
  float k = dot(r,rcentre.xyz);
  // Just antialias the central part of the image.
  float aa0 = float(k > 0.96 ? 3 : k > 0.9 ? 2 : 1);
  float aa = aa0;
  if (iMouse.z > 0.0) aa = 1.0; // Just for comparison
  for (float i = 0.0; i < aa; i++) {
    for (float j = 0.0; j < aa; j++) {
      vec3 col1;
      // What to do with partially visible pixels?
      if (!scene(p,normalize(r+(i-0.5*aa)/aa*drdx+(j-0.5*aa)/aa*drdy),col1)) {
        fragColor = vec4(0,0,0,1);
        return;
      }
      color += col1;
    }
  }
  color /= aa*aa;
  //if (dot(r,rcentre.xyz) > 0.999) color.b = 1.0;
  color *= sqrt(aa0/3.0); // Show aa bands
  if (alert) color.x = 1.0;
  fragColor = vec4(sqrt(color),1);
}

void mainVR(out vec4 fragColor, vec4 eye, vec4 screenpos) {
  eye /= eye.w;
  screenpos /= screenpos.w;
  vec4 rcentre = inverse(iView)*vec4(0,0,-1,0);
  assert(eq(rcentre.w,0.0));
  // The origin of the world is at the viewer, so apply an offset
  // to the camera position (or equivalently, move the scene
  // further away).
  vec3 offset = vec3(0,0,6);
  mainFun(fragColor,eye.xyz+offset,(screenpos-eye).xyz,rcentre.xyz);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
  float scale = 1.0;
  float camera = 4.0;
  vec2 uv = scale*(2.0*fragCoord.xy - iResolution.xy)/iResolution.y;
  vec3 p = vec3(0,0,camera);
  vec3 r = vec3(uv,-2);
  vec3 rcentre = vec3(0,0,-2);
  mainFun(fragColor,p,r,rcentre);
}
