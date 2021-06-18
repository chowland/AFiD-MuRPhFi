### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 54d0751b-c07e-43f7-bb24-ec36793d21a0
push!(LOAD_PATH, pwd())

# ╔═╡ 2a4fd480-6561-11eb-2863-630187ef0d72
using AFiDTools, Plots, HDF5, Statistics, DelimitedFiles, PlutoUI

# ╔═╡ 36f374d0-6561-11eb-12a9-972a3a6d0d33
folder="../test_runs/new_continua"

# ╔═╡ ab177779-9d09-48b3-bf10-38f2d7c8e44c
begin
	fid = h5open(folder*"/outputdir/flowmov/movie_zcut.h5","r")
	lst = keys(fid["temp"])
	varlist = keys(fid)
	close(fid)
end

# ╔═╡ f182d3b0-bd48-11eb-1f79-d79355b36d70
@bind num Slider(parse(Int,lst[1]):parse(Int,lst[end]),show_value=true)

# ╔═╡ 49f6a9ce-a8c5-4777-a8b7-a5d969d71a7d
@bind var Select(varlist)

# ╔═╡ 66a11d40-6561-11eb-1c8a-f1ff8a9d9ab0
T = read_cut(folder,var,num,"z");

# ╔═╡ 7b0ababc-caac-4e7c-97bb-2de977a23e05
# close(fid)

# ╔═╡ 74f21f70-6561-11eb-2297-9fe20d80be5f
grid = read_grid(folder);

# ╔═╡ 9b6df5c0-6561-11eb-06f6-ed301356f103
begin
	gr(dpi=400)
	if size(T)[1]==length(grid.xm)
		heatmap(grid.ym, grid.xm, T, c=:thermal, aspect_ratio=1)
	elseif size(T)[1]==length(grid.xc)
		heatmap(grid.ym, grid.xc, T, c=:balance, aspect_ratio=1)
	else
		heatmap(grid.ymr, grid.xmr, T, c=:ice, aspect_ratio=1)
	end
end

# ╔═╡ Cell order:
# ╟─f182d3b0-bd48-11eb-1f79-d79355b36d70
# ╟─49f6a9ce-a8c5-4777-a8b7-a5d969d71a7d
# ╠═9b6df5c0-6561-11eb-06f6-ed301356f103
# ╠═66a11d40-6561-11eb-1c8a-f1ff8a9d9ab0
# ╠═54d0751b-c07e-43f7-bb24-ec36793d21a0
# ╠═2a4fd480-6561-11eb-2863-630187ef0d72
# ╠═36f374d0-6561-11eb-12a9-972a3a6d0d33
# ╠═ab177779-9d09-48b3-bf10-38f2d7c8e44c
# ╠═7b0ababc-caac-4e7c-97bb-2de977a23e05
# ╠═74f21f70-6561-11eb-2297-9fe20d80be5f
