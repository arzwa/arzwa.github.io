### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ 187a4a1c-4d09-11eb-07ba-1d1f1ba29450
using Plots, StatsPlots, LinearAlgebra, Statistics, Distributions, Measures, 
	QuadGK, Optim

# ╔═╡ 56b80e46-4d44-11eb-1af7-e77b8c42e036
md"""
## Fisher's geometric model
"""

# ╔═╡ 0effa6fc-4d0a-11eb-031a-5b33520807d9
md"""
The distance moved to the optimum upon mutation is 

$$\Delta z = z - z' = z - \sqrt{z^2 - 2rz \cos \theta + r^2}$$

Which Orr states can be approximated by $r \cos \theta - r^2 / 2z$, for which I don't see why... A Taylor series expansion at $r = 0$ gives

$$\approx r \cos \theta - \frac{r^2 \sin^2 \theta }{2z} - \frac{r^3 \sin^2 \theta \cos \theta}{2z^2}$$

$$= r \cos \theta - \frac{r^2}{2z}\Big(\sin^2 \theta - \frac{r \sin^2 \theta \cos \theta}{z}\Big)$$
"""

# ╔═╡ bd46befa-4d08-11eb-39d7-1fd5a959cfd5
Δz(z, θ, r) = z - sqrt(z^2 - 2r*z*cos(θ) + r^2)

# ╔═╡ 0788760c-4d09-11eb-2805-df2e54bf1990
Δz_approx(z, θ, r) = r*cos(θ) - r^2 / 2z

# ╔═╡ 48c3d184-4d09-11eb-13de-85d55e06a64f
begin
	p1 = plot(0:1., θ->Δz(1., θ, 0.1), 
		xlabel="\$\\theta\$", ylabel="\$\\Delta z\$")
	plot!(p1, 0:1., θ->Δz_approx(1., θ, 0.1))
	p2 = plot(0.05:0.05:0.5, r->Δz(1., 1., r))
	plot!(p2, 0.05:0.05:0.5, r->Δz_approx(1., 1., r), 
		xlabel="\$r\$", ylabel="\$\\Delta z\$")
	plot(p1, p2, legend=false, grid=false, size=(500,200))
end

# ╔═╡ e9236234-4d09-11eb-09b8-b37176d02911
md"""
I find this a bit of a weird approximation to pass by without much comment... The approximation is as expected better for small $r$.
"""

# ╔═╡ 70cbb08c-4d0b-11eb-0b99-572a0358a46a
md"""
The crucial change of variables is

$$y = \sqrt{n} \cos \theta$$

bundling the angle and the dimensionality of the phenotypic space in one variable. This change of variables is crucial because the isotropic FGM implies that the density for $y$ is given by a *standard Normal distribution* for sufficiently large $n$. This is an important fact, which makes FGM particularly tractable across different $n$.

Note, here $\theta$ is the angle the mutation vector makes with the line connecting the optimum $O$ and the current phenotype $A$. In other words if a random point on the surface of the $n$-ball around $A$ has angle $\theta$, $\sqrt{n} \cos \theta$ has a standard normal distribution. We can test this: 
"""

# ╔═╡ 70c4da5e-4d0c-11eb-3ec9-f17bab030238
# generates a random vector on the surface of an n-ball
function randunithypersphere(n) 
	x = randn(n)
	return x ./ sum(x .^ 2)
end

# ╔═╡ ef4a9526-4d0c-11eb-2b18-6b3073db8aae
begin
	map([2, 5, 10, 50]) do n
		currz = randn(n)
		spl = map(1:100_000) do i 
			r = randunithypersphere(n)
			# u ⋅ v = ||u||*||v||*cosθ
			# -currz, because we want the angle between currz and r, not O and r
			cosθ = dot(r, -currz)/(norm(-currz) * norm(r))
			√n * cosθ
		end
		plot(Normal(), fill=true, fillalpha=0.2, legend=false)
		density!(spl, color=:black, linewidth=2, 
			title="n = $n", title_loc=:left, titlefont=9, 
			xlabel="\$y\$", ylabel="\$f(y)\$")
	end |> x->plot(x..., grid=false, layout=(1,4), size=(700,120), bottom_margin=3mm)
	# savefig("./img/yisnormal.pdf")
end

# ╔═╡ 40f512d6-4d10-11eb-31ea-5bae119616a6
md"""
This shows that indeed for sufficiently large $n$, $y = \sqrt{n} \cos \theta$, where $\theta$ is the angle of a random point on the $n$-ball centered at $z$ with the line connecting $z$ to the optimum, has a standard normal distribution.
"""

# ╔═╡ b2ee0eba-4d10-11eb-24f5-8dc7aa6e492f
md"""
Now taking Orr's approximate form for $\Delta z$, we obtain on the $y$ scale

$$\Delta z(y) \approx \frac{r}{\sqrt{n}}\Big(y - \frac{r\sqrt{n}}{2z}\Big)$$

as the distance moved towards the optimum upon a mutation of effect $r$ with angle $\theta$ in $n$-dimensional phenotype space. 

*Note again that as a function of $\theta$, the distance moved to the optimum is not dependent on the dimensionality of the phenotype space*! However, we cannot think of $\theta$ easily when $n>2$, so we introduce $y$, which does depend on $n$. While we cannot obtain a density for $\theta$, we can for $y$, which incorporates the dimensionality.

It is now natural to introduce a scaled mutation size 

$$x = \frac{r \sqrt{n}}{2z}$$

So that $\Delta z(y) \approx \frac{2z}{n} x (y - x)$, and a mutation is adaptive as soon as $y > x$. From this we can obtain the probability of a mutation being adaptive easily as

$$\Pr\{\Delta z(y) > 0\} \approx \int_x^\infty f(y) dy = 1 - \Phi(x)$$

**Note** that here the approximation of $\Delta z$ is important, because it allows us to find approximately the threshold value of $y$ for which a mutation is neutral. The expression above is therefore approximate in the sense that the lower limit of the integration is approximate.

"""

# ╔═╡ 1cec3df0-5108-11eb-2dc6-2d3768c60114
begin
	plot(0:0.1:5, Normal(), color=:black, grid=false, size=(350,250), legend=false)
	plot!(0.7:0.1:5, Normal(), fill=true, color=:gray, alpha=0.2,
		xlabel="\$y\$", ylabel="\$f(y)\$")
	annotate!(0.7, 0.02, text("\$x\$", :right))
	annotate!(1.1, 0.07, text("\$P_a\$", :left, :gray))
	# savefig("img/fisherpa1.pdf")
end

# ╔═╡ 9010d1c4-4d3e-11eb-1b32-8d205fd991d7
Δz_y_approx(y, z, n, x) = 2z*x*(y-x)/n

# ╔═╡ ba66891e-4d3e-11eb-1773-2b4e8b562688
Δz_y(y, z, n, x) = z - sqrt(z^2 - 4z^2*x*y/n + 4z^2*x^2/n)

# ╔═╡ 403175d6-4d3f-11eb-2377-63671d52ad70
map([0.1, 0.2, 0.5, 0.8]) do x
	plot(y->Δz_y_approx(y, 1, 10, x), framestyle=:origin)
	plot!(y->Δz_y(y, 1, 10, x), legend=false, grid=false)
	xlabel!("\$y\$")
	ylabel!("\$\\Delta z(y)\$")
	title!("x = $x", title_loc=:left, titlefont=8)
end |> xs -> plot(xs...)

# ╔═╡ 736d5f04-4dea-11eb-2099-7bfa3f2fceb7
begin
	begin
		pp = plot(grid=false, size=(500,300), legend=:outertopright)
		for n=5:5:50
			plot!(pp, 0.01:0.01:1, r->1-cdf(Normal(), r * √n /2), 
				label="n = $n", fg_legend=false)  # z = 1
		end
		plot(pp, xlabel="\$r\$", ylabel="\$P_a\$", 
			fontfamily=:helvetica, guidefont=14, tickfont=9, legendfontsize=10)
	end
	savefig("img/fisherpa.svg")
end

# ╔═╡ 4123c070-4d44-11eb-03e9-5fe04ef8e08f
md"""
This is the first manifestation of the 'cost of complexity', in higher dimensional phenotypic space, for a given mutational effect size, the probability of adaptation decreases markedly.

Making sense of the approximations remains a bit hard, but I should accept them for now and move on...
"""

# ╔═╡ 6d8d654e-4d44-11eb-15e0-e3c3566bd043
md"""
## Orr's model

### The rate of adaptation under Gaussian selection

The expected distance to the optimum upon *adaptive substitution* (i.e. a fixed favorable mutation) is

$$E[z'] = z - E[\Delta z] = z - E[\Delta z | \text{fix}] P [\text{fix}|\text{ad. mutation}] P [\text{ad. mutation}]$$

$$=z - (N \mu \mathrm{d}t) \Pi P_a E[\Delta z| \text{fix}]$$

Where $\Delta z$ is the movement *towards* the optimum. The rate of change of the expected distance to the optimum is

$$\frac{d E[z]}{dt} = -N \mu \Pi P_a E[\Delta z| \text{fix}]$$

**Note:** this ignores potentially deleterious mutations that may fix and increase the distance to the optimum. We are considering the rate of change in the distance to the optimum *during* adaptation.
"""

# ╔═╡ 43b0ee54-4dca-11eb-220d-8f1b1d7dc42d
md"""
The probability of fixing an adaptive mutation is $\Pi \approx 2s$ (this is a result from diffusion theory, see e.g. Rice eq 5.72 for diploid case). Now note that the selection coefficient is defined as

$$\frac{1 + s}{1} = \frac{w(z - Δz_{fav})}{w(z)}$$

Where $\Delta z_{fav}$ is the distance moved to the optimum upon favorable mutation. Now Orr considers Gaussian stabilizing selection, i.e. $w(z) = \exp(-z^2/2)$ so that

$$s = \frac{\exp(-(z - \Delta z_{fav})^2/2)}{\exp(-z^2/2)} - 1 = 
	\exp\Big(\frac{z^2 - (z - \Delta z_{fav})^2}{2}\Big) - 1$$

(note that there is a typo in Orr). For $\Delta z_{fav}$ smallish we have

$$s \approx \exp(z \Delta z_{fav}) - 1 \approx z \Delta z_{fav}$$

So that $\Pi = 2z E[\Delta z|\text{fav}]$ and

$$\frac{d E[z]}{dt} = -2 N \mu P_a z E[\Delta z|\text{fix}] E[\Delta z|\text{fav}]$$

We can measure time on the scale of number of mutations $N \mu$ so that

$$\frac{d E[z]}{dt} = - 2 P_a z E[\Delta z|\text{fix}] E[\Delta z|\text{fav}]$$

Now the goal is to put this in *purely phenotypic terms*, **this is where FGM** comes in, we will derive $E[\Delta z|\text{fix}]$, $E[\Delta z|\text{fav}]$ and $P_a$ from the properties of FGM.
"""

# ╔═╡ 54ff1dc4-4dcb-11eb-0ee5-396b81f325de
md"""
## Integrating FGM

We already know that for sufficiently large $n$, $P_a \approx 1-\Phi(x)$ where $x = \frac{r \sqrt{n}}{2z}$ is Fisher's *standardized unit for the size of a mutation* (which is not that easy to grasp). 

From FGM, we can find that

$$E[\Delta z|\mathrm{fav}] = \frac{r}{\sqrt{n}}
	\frac{\int_x^\infty (y - x) \exp(-y^2/2) dy}{\int_x^\infty \exp(-y^2/2)dy}$$

and

$$E[\Delta z|\mathrm{fix}] = \frac{r}{\sqrt{n}}
	\frac{\int_x^\infty (y - x)^2 \exp(-y^2/2) dy}
         {\int_x^\infty (y - x)\exp(-y^2/2)dy}$$

So that 

$$E[\Delta z|\text{fav}] E[\Delta z|\text{fix}] = \frac{r^2}{n}
	\frac{\int_x^\infty (y - x)^2 \exp(-y^2/2) dy}
         {\int_x^\infty \exp(-y^2/2)dy}$$

and recalling that $P_a \approx 1-\mathrm{\Phi}(x)$

$$\frac{dE[z]}{dt} = -2 P_a z E[\Delta z|\text{fav}] E[\Delta z|\text{fix}] \approx
	- \frac{2z r^2}{n\sqrt{2\pi}}\int_x^\infty (y - x)^2 \exp(-y^2/2) dy$$

Which Orr writes as

$$\frac{dE[z]}{dt} = \frac{2r^2}{n} M z$$

(on the timescale of $(N\mu)^{-1}$ that is)
"""

# ╔═╡ 6654d53a-5022-11eb-3174-a5648a51f4cd
function M(x; upper=1000) 
	integral, err = quadgk(y->(y - x)^2 * exp(-y^2/2), x, 1000)
	integral / √(2π)
end

# ╔═╡ 9db49e14-4de1-11eb-2bcb-1738238b80ad
function dzdt(r, n, z; upper=1000) 
	x = r * √n / (2 * z)
	integral = M(x, upper=upper)
	dzdt = -2z*r^2 * integral / n
end

# ╔═╡ 59b5e692-4de2-11eb-0e1a-bbdf720d667f
begin
	p = plot()
	for n=5:5:50
		plot!(0.01:0.01:1.1, r->-dzdt(r, n, 1.), color=:gray)
		annotate!(0.92, -dzdt(0.9, n, 1.), text("n = $n", :left, 9))
		xmax = optimize(r->dzdt(r, n, 1.), 0., 1.)
		scatter!([xmax.minimizer], [-xmax.minimum], color=:firebrick)
	end
	plot(plot(p, xlabel="\$r\$", ylabel="\$-dz/dt\$",
		legend=false, yscale=:log10, grid=false, size=(300,250)),
	plot(0:0.01:5, x->M(x), grid=false, size=(600,200), legend=false, 
		color=:black, xlabel="\$x\$", ylabel="\$M(x)\$"))
	# savefig("img/dzdt.pdf")
end

# ╔═╡ 238b134e-4de7-11eb-1040-1f4fa20902e0
md"""
This shows Fisher's argument that in higher dimensions, the optimal mutational effect size becomes smaller. Also important is that adaptation is slower globally the higher the dimension of phenotype space. 
"""

# ╔═╡ de887dee-4e0a-11eb-2140-25fad0c58b79
md"""
Orr moves from the dynamics of phenotype $z$ to the dynamics of fitness $w$. Of course, 

$$\frac{\partial E[w(z)]}{\partial t} = \frac{\partial E[w(z)]}{\partial E[z(t)]} \ \frac{\partial E[z(t)]}{\partial t}$$

With Gaussian fitness

$$\frac{\partial E[w(z)]}{\partial E[z(t)]} = -E[z(t)] \exp(-E[z]^2/2)$$

So, writing $w$ for $E[w]$ and $z$ for $E[z]$,

$$\frac{\partial w(z)}{\partial t} = z \exp(-z^2/2) \frac{2 z r^2}{n\sqrt{2\pi}}\int_x^\infty (y - x)^2 \exp(-y^2/2)$$ 
"""

# ╔═╡ 5e921efe-5026-11eb-15c1-f5e730178160
function dwdt(r, n, w; upper=1000)
	x = r * √n / (2 * √(-2 * log(w)))
	integral = M(x, upper=upper)
	dwdt = -4 * w * log(w) * r^2 * integral / n
end

# ╔═╡ a47a100e-5026-11eb-1a99-d1c46820ed3e
let
	p1 = plot(yscale=:log10, xlabel="\$r\$", ylabel="\$-dz/dt\$", 
		title="Rate of phenotypic change", titlefont=10)
	for n=10:5:50
		plot!(0.01:0.01:1.1, r->-dzdt(r, n, 1.), color=:gray)
		n < 30 && annotate!(0.92, -dzdt(0.9, n, 1.), text("n = $n", :left, 9))
		xmax = optimize(r->dzdt(r, n, 1.), 0., 1.)
		scatter!([xmax.minimizer], [-xmax.minimum], color=:firebrick)
	end
	p2 = plot(xlabel="\$r\$", ylabel="\$d\\bar{w}/dt\$",
		title="Rate of fitness increase", titlefont=10)
	for n=10:5:50
		plot!(0.01:0.01:5.1, r->dwdt(r, n, .1), color=:gray)
		xmax = optimize(r->-dwdt(r, n, .1), 0., 5.)
		scatter!([xmax.minimizer], [-xmax.minimum], color=:firebrick)
		n < 30 && annotate!(xmax.minimizer+0.3, 
			-xmax.minimum, text("n = $n", :left, 9))
	end
	p3 = plot(0:0.001:3, x->M(x), grid=false, size=(600,200), legend=false, 
		color=:black, xlabel="\$x\$", ylabel="\$M(x)\$", title="\$M(x)\$")
	annotate!(p3, 1.5, 0.3, text("\$x \\propto r\$", :left))
	annotate!(p3, 1.5, 0.25, text("\$x \\propto \\sqrt{n}\$", :left))
	plot(p3, p1, p2, layout=(1,3), grid=false, legend=false, 
		size=(750,200), fontfamily=:helvetica)
	# savefig("img/dzdt.svg")
end

# ╔═╡ a02e67d0-5183-11eb-27a0-47fbda68b874
plot

# ╔═╡ d042e30a-504e-11eb-37d9-cfffb56876c6
function dwxdt(x, n, w; upper=1000)
	integral = M(x, upper=upper)
	r = x*(2 * √(-2 * log(w))) / √n
	dwdt = -4 * w * log(w) * r^2 * integral / n
end

# ╔═╡ dd83ca16-504e-11eb-2c9f-3bb04603627b
let
	p = plot(xlabel="\$x\$", ylabel="\$d\\bar{w}/dt\$", left_margin=2mm)
	for n=10:5:50
		plot!(p, 0:0.01:5, x->dwxdt(x, n, 0.1), label="n = $n")
	end
	vline!([0.94], label="", color=:black, linestyle=:dash)
	plot(p, legend=:outertopright, size=(350,200), grid=false, fg_legend=:transparent)
	# savefig("img/dwxdt.pdf")
end

# ╔═╡ Cell order:
# ╠═187a4a1c-4d09-11eb-07ba-1d1f1ba29450
# ╟─56b80e46-4d44-11eb-1af7-e77b8c42e036
# ╟─0effa6fc-4d0a-11eb-031a-5b33520807d9
# ╠═bd46befa-4d08-11eb-39d7-1fd5a959cfd5
# ╠═0788760c-4d09-11eb-2805-df2e54bf1990
# ╠═48c3d184-4d09-11eb-13de-85d55e06a64f
# ╟─e9236234-4d09-11eb-09b8-b37176d02911
# ╟─70cbb08c-4d0b-11eb-0b99-572a0358a46a
# ╠═70c4da5e-4d0c-11eb-3ec9-f17bab030238
# ╠═ef4a9526-4d0c-11eb-2b18-6b3073db8aae
# ╟─40f512d6-4d10-11eb-31ea-5bae119616a6
# ╠═b2ee0eba-4d10-11eb-24f5-8dc7aa6e492f
# ╠═1cec3df0-5108-11eb-2dc6-2d3768c60114
# ╠═9010d1c4-4d3e-11eb-1b32-8d205fd991d7
# ╠═ba66891e-4d3e-11eb-1773-2b4e8b562688
# ╠═403175d6-4d3f-11eb-2377-63671d52ad70
# ╠═736d5f04-4dea-11eb-2099-7bfa3f2fceb7
# ╟─4123c070-4d44-11eb-03e9-5fe04ef8e08f
# ╟─6d8d654e-4d44-11eb-15e0-e3c3566bd043
# ╠═43b0ee54-4dca-11eb-220d-8f1b1d7dc42d
# ╠═54ff1dc4-4dcb-11eb-0ee5-396b81f325de
# ╠═6654d53a-5022-11eb-3174-a5648a51f4cd
# ╠═9db49e14-4de1-11eb-2bcb-1738238b80ad
# ╠═59b5e692-4de2-11eb-0e1a-bbdf720d667f
# ╟─238b134e-4de7-11eb-1040-1f4fa20902e0
# ╠═de887dee-4e0a-11eb-2140-25fad0c58b79
# ╠═5e921efe-5026-11eb-15c1-f5e730178160
# ╠═a47a100e-5026-11eb-1a99-d1c46820ed3e
# ╠═a02e67d0-5183-11eb-27a0-47fbda68b874
# ╠═d042e30a-504e-11eb-37d9-cfffb56876c6
# ╠═dd83ca16-504e-11eb-2c9f-3bb04603627b
