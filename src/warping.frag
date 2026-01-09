precision highp float;                 // 浮動小数点精度を高く設定（演算誤差を低減）

uniform float uTime;                  // 経過時間
uniform vec2  u_resolution;            // 画面解像度（幅, 高さ）
uniform float radius;                  // 半径（最終的に色のゆらぎで使用）
uniform float num;                     // 球面上の点の総数
uniform float an;                      // カメラの回転角
out vec4 fragColor;                    // 出力するピクセルの最終 RGBA 値

// Inigo Quilez - Golden spiral sphere mapping demo

// p: 球面上の法線ベクトル（単位ベクトル）
// 戻り値: x = 最も近い黄金螺旋点のインデックス, y = その点との距離
vec2 inverseSF( vec3 p ) { 
    const float kTau = 6.28318530718;       // 2π
    const float kPhi = (1.0 + sqrt(5.0)) / 2.0; // 黄金比 Φ

    // 層（k）を推定
    // 高さ p.z から黄金螺旋ライン上のどの層にいるかを対数式で近似
    float k  = max(
        2.0, 
        floor(
            log2(num * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z))
            / log2(kPhi + 1.0)
        )
    );

    // フィボナッチ近似を使って F(k), F(k+1) を推定
    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));   // (F(k), F(k+1))

    // 黄金螺旋グリッドの係数を生成
    vec2  ka = 2.0 * F / num;                       // 緯度方向の係数
    vec2  kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0)); // 経度方向の係数

    // 逆行列により点 p に最も近い格子番号 c を推定
    mat2 iB = mat2(ka.y, -ka.x, kb.y, -kb.x)
             / (ka.y * kb.x - ka.x * kb.y);

    vec2 c = floor(iB * vec2(
        atan(p.y, p.x),                  // 経度 φ
        p.z - 1.0 + 1.0 / num            // 緯度方向補正
    ));

    // 近傍の 4 点をチェックして最短距離の点を探す
    float d = 8.0;                       // 最短距離（初期値）
    float j = 0.0;                       // 最も近い点の index

    for (int s = 0; s < 4; s++) {
        vec2  uv = vec2(s & 1, s >> 1);  // (0,0),(1,0),(0,1),(1,1)
        float id = clamp(dot(F, uv + c), 0.0, num - 1.0);

        float phi = kTau * fract(id * kPhi);       // 経度
        float cosTheta = 1.0 - (2.0 * id + 1.0) / num; // 緯度 cosθ
        float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

        vec3 q = vec3(
            cos(phi) * sinTheta,
            sin(phi) * sinTheta,
            cosTheta
        );

        float tmp = dot(q - p, q - p);  // 距離の二乗
        if (tmp < d) {
            d = tmp;
            j = id;
        }
    }

    return vec2(j, sqrt(d));            // 点 index と距離
}

// 0~1 の乱数を返す簡易 hash 関数
float hash1(float n) { 
    return fract(sin(n) * 158.5453123);
}

const mat2 m = mat2( 0.80,  0.60, -0.60,  0.80 );

float noise( in vec2 p )
{return sin(p.x)*sin(p.y);}

float fbm4( vec2 p )
{
    float f = 0.0;
    f += 0.5000*noise( p ); p = m*p*2.02;
    f += 0.2500*noise( p ); p = m*p*2.03;
    f += 0.1250*noise( p ); p = m*p*2.01;
    f += 0.0625*noise( p );
    return f/0.9375;
}

float fbm6( vec2 p )
{
    float f = 0.0;
    f += 0.500000*(0.5+0.5*noise( p )); p = m*p*2.02;
    f += 0.250000*(0.5+0.5*noise( p )); p = m*p*2.03;
    f += 0.125000*(0.5+0.5*noise( p )); p = m*p*2.01;
    f += 0.062500*(0.5+0.5*noise( p )); p = m*p*2.04;
    f += 0.031250*(0.5+0.5*noise( p )); p = m*p*2.01;
    f += 0.015625*(0.5+0.5*noise( p ));
    return f/0.96875;
}

vec2 fbm4_2( vec2 p )
{return vec2(fbm4(p), fbm4(p+vec2(7.8)));}

vec2 fbm6_2( vec2 p )
{return vec2(fbm6(p+vec2(16.8)), fbm6(p+vec2(11.5)));}

float func( vec2 q, out vec4 ron )
{
    q += 0.03 * sin( vec2(0.27,0.23) + length(q)*vec2(4.1,4.3));

    vec2 o = fbm4_2( 0.9*q );
    o += 0.04 * sin( vec2(0.12,0.14) + length(o) );
    vec2 n = fbm6_2( 3.0*o );
    ron = vec4( o, n );

    float f = 0.5 + 0.5*fbm4( 1.8*q + 6.0*n );
    return mix( f, f*f*f*3.5, f*abs(n.x) );
}

//---------------------------------------------------------------
// main
//---------------------------------------------------------------
void main() {

    // ピクセル座標を正規化（中心を (0,0) に）
    vec2 fragCoord = gl_FragCoord.xy;
    vec2 p = (-u_resolution.xy + 2.0 * fragCoord.xy) / u_resolution.y;

    //カメラ設定
    //カメラ位置。anにより、カメラが円周上を回る。高さは1.0
    vec3 ro = vec3(2.5 * cos(an), 1.0, 2.5 * sin(an));
    vec3 ta = vec3(0.0, 1.0, 0.0); // 注視点。カメラはここに向かっている

    vec3 ww = normalize(ta - ro); //カメラの前方向(正規化したta-ro)
    vec3 uu = normalize(cross(ww, vec3(0.0, 1.0, 0.0))); //右方向ベクトル(wwとy軸との外積)
    vec3 vv = normalize(cross(uu, ww));  //上方向ベクトル(uuとwwの外積)
    vec3 rd = normalize(p.x * uu + p.y * vv + 1.5 * ww); 

    // 初期化
    vec3 sc = vec3(0.0, 1.0, 0.0); // 球の中心
    vec3 col = vec3(1.0);
    float tmin = 10000.0;         // 最小の交差距離(レイの最近接交差のt値(レイがどれだけ進んだか))
    //大きい初期値を設定しているのは、レイが何も当たらない場合を示すため
    vec3  nor = vec3(0.0);        //法線
    float occ = 1.0;              //影の係数
    vec3  pos = vec3(0.0);        //交差点

    // 地面 y=0 との交差判定
    float h = (0.0 - ro.y) / rd.y;
    if (h > 0.0) {
        tmin = h; //最初の交差距離
        nor = vec3(0.0, 1.0, 0.0); //地面の法線(上向き)
        pos = ro + h * rd; //交点座標の計算
        vec3 di = sc - pos; //地面上の交点から球中心までのベクトル(球との距離や影処理に使用)
        float l = length(di);
        occ = 1.0 - dot(nor, di / l) / (l * l);
        col = vec3(1.0);
    }

    // 球との交差判定（2次方程式）
    vec3  ce = ro - sc;
    float b = dot(rd, ce);
    float c = dot(ce, ce) - 1.0;   
    h = b * b - c;                 //判別式

    if (h > 0.0) {
    h = -b - sqrt(h);          
    if (h < tmin) {
        tmin = h;
        nor = normalize(ro + h * rd - sc); 
        occ = 0.5 + 0.5 * nor.y;           
    }

    vec2 fi = inverseSF(nor);
    float d = fi.y; // 最近傍点までの距離

    // 球面座標ベースの滑らかなUV
    float phi   = atan(nor.z, nor.x);                     
    float theta = acos(clamp(nor.y, -1.0, 1.0));          
    vec2 uv = vec2(phi / 6.28318, theta / 3.14159);       

    // numに応じて模様の細かさを変える
    const float baseNum = 500.0;
    float detail = sqrt(clamp(num / baseNum, 0.2, 5.0));

    vec2 q = (uv - 0.5) * (3.0 * detail);

    // fbm ドメインワーピングで模様
    vec4 ron;
    float f = func(q, ron);

    vec3 flowCol = vec3(0.0);
    flowCol = mix( vec3(0.15,0.10,0.25), vec3(0.35,0.10,0.08), f );
    flowCol = mix( flowCol, vec3(0.9,0.9,0.9), clamp(dot(ron.zw, ron.zw), 0.0, 1.0) );

    // fi.y を使って構造のエッジを強調
    float edge = exp(-d * (40.0 + 0.03 * num)); 
    // フィボナッチ点の近くが強くなる

    // fbmのコントラストを強める
    float sharpF = smoothstep(0.2, 0.8, f);

    // flowColに追加
    flowCol *= 0.8 + 0.4 * sharpF;     // fbm のコントラスト
    flowCol *= 0.9 + 0.3 * edge;       // フィボナッチ構造の強調

    vec3 baseCol = flowCol;

    // フィボナッチ点の存在感を少し足す
    /*
    float halo = 1.0 - smoothstep(0.0, 0.06, d); // 点の近くほど 1 に近い
    float contrast = mix(1.0, 1.25, halo);
    */

    col = baseCol;

    // 半径と距離で揺らぎ
    float wave = 1.0 + 0.04 * sin(60.0 * radius * d);
    col *= wave;

    col *= 1.2;
}

    // 背景フェード
    if (tmin < 100.0) {
        pos = ro + tmin * rd;
        col *= occ;
        col = mix(col, vec3(1.0), 1.0 - exp(-0.003 * tmin * tmin));
    }

    // ガンマ補正
    col = sqrt(col);
    fragColor = vec4(col, 1.0);
}
