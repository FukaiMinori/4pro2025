precision highp float;

out vec4 fragColor;
uniform vec2 u_resolution;
uniform float radius;
uniform float alpha;

void main(){
    vec2 p = gl_FragCoord.xy/u_resolution.x;

    if(length(p)<radius){
        fragColor = vec4(1.0,1.0,0.0,alpha);
    }
}