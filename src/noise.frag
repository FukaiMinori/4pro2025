precision highp float;

out vec4 fragColor;

uniform vec2 u_resolution;
uniform float power;

const uint UINT_MAX = 0xffffffffu;
uvec3 k = uvec3(0x456789abu, 0x6789ab45u, 0x89ab4567u);
uvec3 u = uvec3(1, 2, 3);
uvec2 uhash22(uvec2 n){
    n ^= (n.yx << u.xy);
    n ^= (n.yx >> u.xy);
    n *= k.xy;
    n ^= (n.yx << u.xy);
    return n * k.xy;
}
uvec3 uhash33(uvec3 n){
    n ^= (n.yzx << u);
    n ^= (n.yzx >> u);
    n *= k;
    n ^= (n.yzx << u);
    return n * k;
}
vec2 hash22(vec2 p){
    uvec2 n = floatBitsToUint(p);
    return vec2(uhash22(n)) / vec2(UINT_MAX);
}
vec3 hash33(vec3 p){
    uvec3 n = floatBitsToUint(p);
    return vec3(uhash33(n)) / vec3(UINT_MAX);
}
float hash21(vec2 p){ //二次元ベクトルに対してフロートを一つ返す
    uvec2 n = floatBitsToUint(p);
    return float(uhash22(n).x) / float(UINT_MAX);
    //nesting approach
    //return float(uhash11(n.x+uhash11(n.y)) / float(UINT_MAX)
}
float hash31(vec3 p){
    uvec3 n = floatBitsToUint(p);
    return float(uhash33(n).x) / float(UINT_MAX);
    //nesting approach
    //return float(uhash11(n.x+uhash11(n.y+uhash11(n.z))) / float(UINT_MAX)
}

float h3(float x){
    return ((3.0 -2.0*x)*x*x);
}

vec4 hokan(vec2 p){
    vec4 bl = vec4(vec3(hash21(floor(p))), 1.0);
    vec4 br = vec4(vec3(hash21(floor(p) + vec2(1,0))), 1.0);
    vec4 tl = vec4(vec3(hash21(floor(p) + vec2(0,1))), 1.0);
    vec4 tr = vec4(vec3(hash21(floor(p) + vec2(1,1))), 1.0);

    vec4 bot = mix(bl, br, fract(p.x));
    vec4 top = mix(tl, tr, fract(p.x));
    vec4 v = mix(bot, top, fract(p.y));
    return v;
}

float vnoise(vec2 p){
    float bl = hash21(floor(p));
    float br = hash21(floor(p) + vec2(1,0));
    float tl = hash21(floor(p) + vec2(0,1));
    float tr = hash21(floor(p) + vec2(1,1));

    float bot = mix(bl, br, fract(p.x));
    float top = mix(tl, tr, fract(p.x));
    float v = mix(bot, top, fract(p.y));
    return v;
}

float fbm(vec2 p){
    float val = vnoise(p);
    float mag = 2.0;
    float g = 2.0;

    for(int i=1; i<20; i++){
        val += vnoise(mag*p)/g;
        mag *= 2.0;
        g *= 2.0;
    }
    return 0.5+0.5*val;
}

float warp(vec2 p){
    vec2 q = vec2(fbm(p),fbm(p + vec2(5.2, 1.3)));
    vec2 r = vec2(fbm(p+ 2.0*q + vec2(1.7, 9.2)),fbm(p+ 2.0*q + vec2(8.3, 2.8)));

    return vnoise(p + 4.0*r);
}


void main(){
    vec2 p = gl_FragCoord.xy / u_resolution.x;
    p = 20.0*p;
    int ipx = int(floor(p.x));
    int ipy = int(floor(p.y));
    //p = fract(p);

    fragColor = vec4(vec3(0.5+power*warp(p)),1.0);

    //fragColor = vec4(vec3(hash21(floor(p))), 1.0);

    /*
    vec4 colLB = vec4(0.0, 0.3059, 0.9216, 1.0);
    vec4 colRB = vec4(0.9843, 0.6078, 0.1176, 1.0);
    vec4 colLT = vec4(0.1098, 0.9098, 1.0, 1.0);
    vec4 colRT = vec4(0.4431, 1.0, 0.502, 1.0);

    vec4 colB;
    vec4 colT;

    if(ipx%2==0 && ipy%2==0){
        colB = mix(colLB, colRB, h3(p.x));
        colT = mix(colLT, colRT, h3(p.x));
        fragColor = mix(colB, colT, h3(p.y));
    }
    else if(ipx%2==1 && ipy%2==0){
        colB = mix(colRB, colLB, h3(p.x));
        colT = mix(colRT, colLT, h3(p.x));
        fragColor = mix(colB, colT, h3(p.y));
    }
    else if(ipx%2==0 && ipy%2==1){
        colT = mix(colLB, colRB, h3(p.x));
        colB = mix(colLT, colRT, h3(p.x));
        fragColor = mix(colB, colT, h3(p.y));
    }
    else if(ipx%2==1 && ipy%2==1){
        colT = mix(colRB, colLB, h3(p.x));
        colB = mix(colRT, colLT, h3(p.x));
        fragColor = mix(colB, colT, h3(p.y));
    }
    */
    
}