precision highp float;

uniform vec2 u_resolution;
uniform float theta;
uniform float uTime;

out vec4 fragColor;



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

float mindist(vec2 p){
    vec2 f = floor(p);
    float m = 10.0;
    for(int i=-5; i<=5; i++){
        for(int j=-5; j<=5; j++){
            vec2 l = f+vec2(i,j);
            m = min(m, distance(l+theta*hash22(l),p));
        }
    }
    return m;
}

float cellular(vec2 p){
    vec2 f = floor(p);
    float m1 = 10.0;
    float m2 = 10.0;
    float m3 = 10.0;
    float m4 = 10.0;
    for(int i=-5; i<=5; i++){
        for(int j=-5; j<=5; j++){
            vec2 l = f+vec2(i,j);
            float d = distance(l+sin(uTime)+theta*hash22(l),p);
            if(d<m1){
                m4 = m3;
                m3 = m2;
                m2 = m1;
                m1 = d;
            }else if(d<m2){
                m4 = m3;
                m3 = m2;
                m2 = d;
            }else if(d<m3){
                m4 = m3;
                m3 = d;
            }else if(d<m4){
                m4 = d;
            }
        }
    }
    return 0.3*m1-0.2*m2-0.2*m3+0.5*m4;
}

void main(){
    vec2 p = gl_FragCoord.xy / u_resolution.x;

    p = 5.0*p;
    fragColor = vec4(vec3(cellular(p)),1.0);
}