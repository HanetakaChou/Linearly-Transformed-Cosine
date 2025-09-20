void main(
	in uint d3d_VertexID: SV_VertexID,
	out float4 d3d_Position : SV_POSITION,
	out float2 out_uv : TEXCOORD0
)
{
	const float2 positions[3] = { float2(-1.0, -3.0), float2(-1.0, 1.0), float2(3.0, 1.0) };
	const float2 uv[3] = { float2(0.0, 2.0), float2(0.0, 0.0), float2(2.0, 0.0) };

	out_uv = uv[d3d_VertexID];
	d3d_Position = float4(positions[d3d_VertexID], 0.5, 1.0);
}