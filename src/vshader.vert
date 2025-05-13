precision mediump float;

in vec3 vPosition;

void main(void) {
  //gl_Position = vec4(vPosition, 1.0);
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}
