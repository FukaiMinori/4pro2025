//import './App.css'
import * as THREE from 'three'
import { Canvas, extend } from '@react-three/fiber'
import { Hud, OrthographicCamera, shaderMaterial } from '@react-three/drei'
import { useControls } from 'leva'
//import './styles.css'
import vertexShader from "./vshader.vert?raw";
import shaderA from "./sphericalFibo9.frag?raw";
import shaderB from "./sphericalFibo11.frag?raw";
import shaderC from "./sphericalFibo12.frag?raw";
import { useEffect, useState } from 'react';
import { Box, Select, MenuItem } from "@mui/material";

declare global
{ namespace JSX
 { interface IntrinsicElements
        { "diffMaterial": any
        }
 }
}

function App() {
  const [shaderKey, setShaderKey] = useState<any>("c");
  const { num, an, embedSmall, embedSmall2} = useControls({
    num :{value:130.0, min: 10.0, max: 1000.0, step:1.0},
    an :{value:0.0, min:-3.14159 / 2.0, max:3.14159 / 2.0, step:0.1},
    embedSmall :{value:1.7, min: 0.5, max: 2.0, step:0.005},
    embedSmall2 :{value:1.7, min: 0.5, max: 2.0, step:0.005}
  });

  const [uTime, setUTime] = useState(0.0);

  const shaders:any = {
    a:{
      name:"recursive",
      fragment:shaderA,
    },
    b:{
      name:"celluler noise",
      fragment:shaderB,
    },
    c:{
      name:"voronoi",
      fragment:shaderC,
    },
  };

  useEffect(()=>{
    const interval = setInterval(()=>{
      setUTime((uTime)=>uTime+0.1);
      //console.log(uTime);
    },100)
    return ()=>clearInterval(interval);
  },[uTime]
    )
  return (
    <>
      <Box sx={{ width: "100vw", height: "100vh" }}>
    {/* UI */}
    <Box sx={{ position: "absolute", zIndex: 10, p: 2 }}>
      <Select
        value={shaderKey}
        onChange={(e) => {console.log(e.target.value);setShaderKey(e.target.value)}}
      >
        {Object.entries(shaders).map(([key, s]:any) => (
          <MenuItem key={key} value={key}>
    {s.name}
          </MenuItem>
        ))}
      </Select>
  </Box>
	<div style={{ height: "100dvh", width: "100dvw" }}>
    <Scene num={num} an={an} uTime={uTime} embedSmall={embedSmall} embedSmall2={embedSmall2} shaderKey={shaderKey} shader={shaders[shaderKey].fragment}/>
	</div>
    </Box>
  </>
  )
}

function Scene({num, an, uTime, embedSmall, embedSmall2, shader, shaderKey}:{num:number, an:number, uTime:number, embedSmall:number, embedSmall2:number, shader:any, shaderKey:any}){
  return(
    <Canvas>
      <Hud>
	<OrthographicCamera
	  makeDefault
	  top={1}
	      right={1}
	      bottom={-1}
	      left={-1}
	      near={0}
	      far={1}
	/>
	<OrthographicCamera  makeDefault top={1} right={1} bottom={-1} left={-1} near={0} far={1}/>
	<mesh>
	  <planeGeometry args={[2,2]}/>
	  {/* @ts-ignore TS2339: Property 'diffMaterial' does not exist on type 'JSX.IntrinsicElements'.*/}
	  <diffMaterial key={shaderKey} fragmentShader={shader} glslVersion={THREE.GLSL3} num={num} an={an} uTime={uTime} embedSmall={embedSmall} embedSmall2={embedSmall2}/>
	</mesh>
      </Hud>
    </Canvas>
  )
}
const DiffMaterial = shaderMaterial(
  {u_resolution: new THREE.Vector2(window.innerWidth, window.innerHeight),
    num: 130.0,
    an: 0.0,
    embedSmall: 1.7,
    embedSmall2: 1.7,
    uTime: 0.0
   },
  vertexShader,
  shaderC
);
extend({DiffMaterial});


export default App
