SamplerState clamp_sampler : register(s0);

Texture2D in_color : register(t0);

float3 aces_fitted(float3 color);

float3 ToSRGB(float3 v);

void main(
	in float4 d3d_Position : SV_POSITION,
	in float2 in_uv : TEXCOORD0,
	out float4 out_color : SV_TARGET0
	)
{
	float3 col = in_color.Sample(clamp_sampler, in_uv).rgb;

	col = aces_fitted(col);

	col = ToSRGB(col);

	out_color = float4(col, 1.0);
}

float3 rrt_odt_fit(float3 v)
{
	float3 a = v * (v + 0.0245786) - 0.000090537;
	float3 b = v * (0.983729 * v + 0.4329510) + 0.238081;
	return a / b;
}

float3 aces_fitted(float3 color)
{
	float3x3 ACES_INPUT_MAT = 
		float3x3(
			float3(0.59719, 0.35458, 0.04823), // row 0
			float3(0.07600, 0.90834, 0.01566), // row 1
			float3(0.02840, 0.13383, 0.83777)  // row 2
			);

	float3x3 ACES_OUTPUT_MAT = 
		float3x3(
			float3(1.60475, -0.53108, -0.07367), // row 0
			float3(-0.10208, 1.10813, -0.00605), // row 1
			float3(-0.00327, -0.07276, 1.07602)  // row 2
			);

	color = mul(ACES_INPUT_MAT, color);

	// Apply RRT and ODT
	color = rrt_odt_fit(color);

	color = mul(ACES_OUTPUT_MAT, color);

	// Clamp to [0, 1]
	color = saturate(color);

	return color;
}

float3 ToSRGB(float3 v)
{
	return pow(v, (1.0 / 2.2));
}