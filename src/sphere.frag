precision highp float;

uniform float theta;
uniform float N;
uniform float uTime;
uniform vec2  u_resolution;   // Shadertoy互換名に統一
out vec4 fragColor;

// ---------------------------------------------------------------
// Inigo Quilez - https://iquilezles.org/
// Golden spiral sphere mapping demo
// ---------------------------------------------------------------

vec2 inverseSF( vec3 p , float kNum) {
    const float kTau = 6.28318530718;
    const float kPhi = (1.0 + sqrt(5.0)) / 2.0;
    //const float kNum = N;

    float k  = max(2.0, floor(log2(kNum * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z)) / log2(kPhi + 1.0)));
    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));

    vec2  ka = 2.0 * F / kNum;
    vec2  kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0));

    mat2 iB = mat2(ka.y, -ka.x, kb.y, -kb.x) / (ka.y * kb.x - ka.x * kb.y);
    vec2 c = floor(iB * vec2(atan(p.y, p.x), p.z - 1.0 + 1.0 / kNum));

    float d = 8.0;
    float j = 0.0;
    for (int s = 0; s < 4; s++) {
        vec2  uv = vec2(s & 1, s >> 1);
        float id = clamp(dot(F, uv + c), 0.0, kNum - 1.0);

        float phi = kTau * fract(id * kPhi);
        float cosTheta = 1.0 - (2.0 * id + 1.0) / kNum;
        float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

        vec3 q = vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
        float tmp = dot(q - p, q - p);
        if (tmp < d) {
            d = tmp;
            j = id;
        }
    }
    return vec2(j, sqrt(d));
}


vec4 nearestSFids( vec3 p , float kNum) {
    const float kTau = 6.28318530718;
    const float kPhi = (1.0 + sqrt(5.0)) / 2.0;

    float k  = max(2.0, floor(log2(kNum * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z)) / log2(kPhi + 1.0)));
    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));

    vec2  ka = 2.0 * F / kNum;
    vec2  kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0));

    mat2 iB = mat2(ka.y, -ka.x, kb.y, -kb.x) / (ka.y * kb.x - ka.x * kb.y);
    vec2 c = floor(iB * vec2(atan(p.y, p.x), p.z - 1.0 + 1.0 / kNum));

    float d = 8.0;
    float j = 0.0;
    vec3 pt = p;
    vec3 q;
    vec4 ids = vec4(0);
    for (int s = 0; s < 4; s++) {
        vec2  uv = vec2(s & 1, s >> 1);
        float id = clamp(dot(F, uv + c), 0.0, kNum - 1.0);

        float phi = kTau * fract(id * kPhi);
        float cosTheta = 1.0 - (2.0 * id + 1.0) / kNum;
        float sinTheta = sqrt(1.0 - cosTheta * cosTheta);
	ids[s] = id;
    }
    return ids;
}

float hash1(float n) { return fract(sin(n) * 158.5453123); }


vec3 sphericalPt(vec3 center, float rad, vec3 camera, vec3 ray){
  vec3 v = center - camera;
  float b = dot(ray, v);
  float c = dot(v,v) - rad*rad;
  if (b*b - c > 0.0){
    float d = b + sqrt(b*b - c);
    return(camera + d*ray);
  }
  else{
    return(vec3(100000.0));
  }
}

void main() {
    vec2 fragCoord = gl_FragCoord.xy;
    vec2 p = (-u_resolution.xy + 2.0 * gl_FragCoord.xy) / u_resolution.y;
    float an = 0.5 * uTime;

    // camera movement
    //float an = 0.5 * iTime;
    vec3 ro = vec3(5.0 * cos(theta), 1.0, 5.0 * sin(theta));
    vec3 ta = vec3(0.0, 1.0, 0.0);
    vec3 ww = normalize(ta - ro);
    vec3 uu = normalize(cross(ww, vec3(0.0, 1.0, 0.0)));
    vec3 vv = normalize(cross(uu, ww));
    vec3 rd = normalize(p.x * uu + p.y * vv + 1.5 * ww);

    vec3 sc = vec3(0.0, 1.0, 0.0);
    vec3 col = vec3(1.0);

    float tmin = 10000.0;
    vec3  nor = vec3(0.0);
    float occ = 1.0;
    vec3  pos = vec3(0.0);

    // plane
    float h = (0.0 - ro.y) / rd.y;
    if (h > 0.0) {
        tmin = h;
        nor = vec3(0.0, 1.0, 0.0);
        pos = ro + h * rd;
        vec3 di = sc - pos;
        float l = length(di);
        occ = 1.0 - dot(nor, di / l) * 1.0 * 1.0 / (l * l);
        col = vec3(1.0);
    }

    // sphere
    float rad = 2.0;
    vec3 v = sphericalPt(sc, 1.1*rad, ro, rd);
    vec3 np = sc;
    //vec3  ce = ro - sc;
    //float b = dot(rd, ce);
    //float c = dot(ce, ce) - 1.0;
    //h = b * b - c; // discriminant
    //if (h > 0.0) {
    if (length(v) < 10000.0) {
      h = length(v-ro);
      //h = -b - sqrt(h);
      //if (h < tmin) {
      tmin = h;
      nor = normalize(v - sc);
      //occ = 0.5 + 0.5 * nor.y;

      vec2 fi = inverseSF(nor, N);
      np = rad*nearestSF(nor, N) + sc;
      //col *= 1.0 + 0.1 * sin(250.0 * fi.y);
      //col = vec3(1.0,0.0,1.0);
      //col *= 0.0 + 25.0*(fi.y*fi.y);
      //col *= 1.5;
      //col = vec3(1.0,1.0,0.0);
    }

    v = sphericalPt(0.8*np, 0.2, ro, rd);
    if (length(v) < 10000.0) {
      h = length(v-ro);
      if (h < tmin) {
    	tmin = h;
    	nor = normalize(v - np);
    	occ = 0.5 + 0.5 * nor.y;
      }
      col = vec3(1.0,0.0,0.0);
    }
    

    if (length(v) < 10000.0) {
      //if (tmin < 10000.0) {
      //pos = ro + tmin * rd;
      //pos = v;
      col *= occ;
      //col = mix(col, vec3(1.0), 1.0 - exp(-0.003 * tmin * tmin));
    }

    col = sqrt(col);
    fragColor = vec4(col, 1.0);
}
