### A Pluto.jl notebook ###
# v0.15.1

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

# ╔═╡ 5c0baaa4-0d7e-47de-b29e-5bc7281b83cd
using Pkg

# ╔═╡ 4495b957-b8b4-457b-a01a-f3a5d14e4d7c
Pkg.activate()

# ╔═╡ 54d0751b-c07e-43f7-bb24-ec36793d21a0
push!(LOAD_PATH, pwd())

# ╔═╡ 2a4fd480-6561-11eb-2863-630187ef0d72
using AFiDTools, PyPlot, HDF5, Statistics, DelimitedFiles, ColorSchemes, LaTeXStrings, PlutoUI

# ╔═╡ 3c6af843-fd4d-4d8a-a7e0-0e8a18016e3e
128*3

# ╔═╡ 2760ca73-a629-4955-a86e-d4b40718a696
plt.style.use("ggplot")

# ╔═╡ 36f374d0-6561-11eb-12a9-972a3a6d0d33
folder="../examples/SupercooledDiscGrowth"

# ╔═╡ ab177779-9d09-48b3-bf10-38f2d7c8e44c
begin
	fid = h5open(folder*"/outputdir/flowmov/movie_zcut.h5","r")
	lst = keys(fid["temp"])
	varlist = keys(fid)
	close(fid)
end

# ╔═╡ f182d3b0-bd48-11eb-1f79-d79355b36d70
# @bind num Slider(parse(Int,lst[1]):100,show_value=true)
@bind num Slider(parse(Int,lst[1]):parse(Int,lst[end]),show_value=true)

# ╔═╡ 49f6a9ce-a8c5-4777-a8b7-a5d969d71a7d
@bind Var Select(varlist)

# ╔═╡ 66a11d40-6561-11eb-1c8a-f1ff8a9d9ab0
T = read_cut(folder,Var,num,"z");

# ╔═╡ 74f21f70-6561-11eb-2297-9fe20d80be5f
grid = read_grid(folder)

# ╔═╡ 9b6df5c0-6561-11eb-06f6-ed301356f103
begin
	clf()
	Tmax = maximum(abs.(T))
	if size(T)[2]==length(grid.ym)
		pcolormesh(grid.yc, grid.xc, T, cmap="RdBu_r", vmin=-Tmax, vmax=Tmax)
	else
		pcolormesh(grid.ycr, grid.xcr, T, cmap="Blues_r", vmin=0, vmax=1)
	end
	# contour(xm,ym,T,levels=[0.5])
	# contour(ym,xm,T, [0.5])
	colorbar()
	title(Var)
	ax2 = gca()
	ax2.set_aspect("equal")
	gcf()
end

# ╔═╡ d2807e9b-c632-4230-b7c2-69f0923b34e6
begin
	clf()
	if size(T)[1]==length(grid.xmr)
		plot(grid.xmr, T[:,end÷2],"o")
	elseif size(T)[1]==length(grid.xc)
		plot(grid.xc, T[:,end÷2], "o")
	else
		plot(grid.xm, T[:,end÷2], "o")
	end
	gcf()
end

# ╔═╡ 861337de-4bee-4e68-b592-70e1acf2f63d
length(grid.ym), size(T)

# ╔═╡ Cell order:
# ╠═f182d3b0-bd48-11eb-1f79-d79355b36d70
# ╠═49f6a9ce-a8c5-4777-a8b7-a5d969d71a7d
# ╠═9b6df5c0-6561-11eb-06f6-ed301356f103
# ╠═d2807e9b-c632-4230-b7c2-69f0923b34e6
# ╠═3c6af843-fd4d-4d8a-a7e0-0e8a18016e3e
# ╠═861337de-4bee-4e68-b592-70e1acf2f63d
# ╠═66a11d40-6561-11eb-1c8a-f1ff8a9d9ab0
# ╠═5c0baaa4-0d7e-47de-b29e-5bc7281b83cd
# ╠═4495b957-b8b4-457b-a01a-f3a5d14e4d7c
# ╠═54d0751b-c07e-43f7-bb24-ec36793d21a0
# ╠═2a4fd480-6561-11eb-2863-630187ef0d72
# ╠═2760ca73-a629-4955-a86e-d4b40718a696
# ╠═36f374d0-6561-11eb-12a9-972a3a6d0d33
# ╠═ab177779-9d09-48b3-bf10-38f2d7c8e44c
# ╠═74f21f70-6561-11eb-2297-9fe20d80be5f
