### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ fac60950-c64c-4d88-bd27-82a60b737d8f
using PlutoUI, PlutoPlotly, LinearAlgebra, Printf

# ╔═╡ 3373b2f8-d052-4103-96c3-99cbc4d33f43
md"""
# Diseño Jaula de Helmholtz 
### Proyecto Kuntur 

Esta seccion narra el diseño de la bobina cuadrada de helmholtz propuesto por el grupo kuntur para la validacion del control de un CubeSat.

"""

# ╔═╡ 80568fe3-07b0-4144-83a7-ba8700988f33
TableOfContents(title = "Contenido", indent = true)

# ╔═╡ 6f80bf86-6cb4-45be-8cbe-669373fcec5b
md"""
## 0. ¿Qué estamos construyendo y por qué?

 Un **CubeSat** es un satélite pequeño hecho de "cubos" de unos 10 cm. Cada cubo es una **unidad (U)**. Un CubeSat de 1U mide, según el estándar, unos **10 × 10 × 11,35 cm**.

En órbita, el satélite necesita saber hacia dónde apunta y poder girar cuando quiera. A eso se le llama **control de actitud**.Dado que el CubeSats que proponemos estara ubicado en la orbita terrestre baja (LEO) este control se hara apoyándose en el campo magnético de la Tierra el cual:
- se miden con *magnetómetros* (brújulas electrónicas) para saber su orientación,
-  se interactua con *magnetorquers* para girar.

En el laboratorio, donde el campo magnético es fijo, se deben realizar pruebas antes de enviarlo.Para probar sus sensores y actuadores necesitamos fabricar un campo magnético a voluntad: con la magnitud y la dirección que tendría en órbita.

La solución es fabricar una jaula de Helmholtz. La cual se compone de tres pares de bobinas cuadradas, un par por cada eje (X, Y, Z). Cada uno produce un campo uniforme en su dirección, en la zona central. Combinando los tres pares podemos crear un campo en cualquier dirección, y además cancelar el de la Tierra.

"""

# ╔═╡ 41024808-f97d-41f2-9cbe-87fe50712b39
md"""
---
## 1. Homogeneidad: volumen de la zona uniforme

Las bobinas producen un campo muy parejo cerca del centro, que se va deformando a medida que nos alejamos. El satélite completo tiene que quedar dentro de esta zona, sin importar cómo esté girado.

Una forma de garantizarlo es encontrar la esfera más pequeña que contiene al satélite. El diámetro de esa esfera es la diagonal mayor del cubo.

Esa diagonal sale del teorema de Pitágoras en 3D:

```math
D = \sqrt{a^2 + b^2 + c^2}
```

Para el CubeSat de 1U (``a = b = 10`` cm, ``c = 11{,}35`` cm):

```math
D = \sqrt{2\,(10^2) + (11{,}35)^2} = \sqrt{328{,}82} \approx 18{,}13\ \text{cm}
```

"""

# ╔═╡ 060bd874-245b-4d00-9140-f3bc9060b16f
md"""
Largo ``a`` (cm): $(@bind a_cm Slider(5:0.05:40, default = 10, show_value = true))

Ancho ``b`` (cm): $(@bind b_cm Slider(5:0.05:40, default = 10, show_value = true))

Alto ``c`` (cm): $(@bind c_cm Slider(5:0.05:40, default = 11.35, show_value = true))
"""

# ╔═╡ 93b6a733-c764-4525-a68c-1ed211664d6c
# Diagonal mayor de la caja = diámetro de la zona homogénea necesaria (cm)
D_cm = sqrt(a_cm^2 + b_cm^2 + c_cm^2)

# ╔═╡ 8ae398d3-a0cf-46b7-ae5f-807fb9d5fc6d
md"""
---
## 2. Campo necesario: ¿cuántos teslas deben producir las bobinas?

Se parte de dos datos:
- En **Bogotá**, el campo terrestre es de unos **30 µT**.
- En la órbita que se quiere simular (la región del primer cinturón de Van Allen) el campo también es de unos **30 µT**.

Los campos magnéticos son vectores que se suman. Dentro de la jaula, el satélite siente

```math
\vec B_\text{total} = \vec B_\text{Tierra} + \vec B_\text{bobinas}
```

Queremos que ``\vec B_\text{total}`` sea igual al campo de la órbita, ``\vec B_\text{órbita}``. Despejando, las bobinas deben producir

```math
\vec B_\text{bobinas} = \vec B_\text{órbita} - \vec B_\text{Tierra}
```

Las bobinas cancelan el campo de Bogotá y crean el de la órbita. Cuánto campo hace falta depende del ángulo ``\theta`` entre las dos vectores:

```math
|\vec B_\text{bobinas}| = \sqrt{B_\text{Tierra}^2 + B_\text{órbita}^2 - 2\,B_\text{Tierra}\,B_\text{órbita}\cos\theta}
```
"""

# ╔═╡ 4a31f1ed-1a0b-4b71-b1ad-76ed30d9d300
md"""
Campo de la Tierra en el laboratorio (µT): $(@bind B_tierra Slider(20:1:65, default = 30, show_value = true))

Campo de la órbita a simular (µT): $(@bind B_orbita Slider(10:1:65, default = 30, show_value = true))

Ángulo ``\theta`` entre ambos campos (°): $(@bind θ_deg Slider(0:5:180, default = 120, show_value = true))

Factor de seguridad: $(@bind FS Slider(1:0.1:3, default = 1.2, show_value = true))
"""

# ╔═╡ 79557bf3-5dc6-4a08-b56c-fedbff16e618
begin
	# Campo que deben crear las bobinas para el ángulo elegido (µT)
	B_bobinas_θ = sqrt(max(0, B_tierra^2 + B_orbita^2 - 2B_tierra * B_orbita * cosd(θ_deg)))

	# Peor caso: θ = 180° (campos opuestos). Hay que cancelar uno entero y crear el otro.
	B_peor = B_tierra + B_orbita

	# Campo objetivo de diseño, con factor de seguridad (µT, por eje)
	B_obj = FS * B_peor
end

# ╔═╡ 2654b6d3-13c5-47af-a5f8-689ca33ba2c4
md"""
---
## 3. Geometría: tamaño y separación de las bobinas

### 3.1 ¿Cómo crea campo una bobina?

Cuando la corriente recorre una espira, crea un campo magnético que la atraviesa. Su dirección sigue la regla de la mano derecha: los dedos en direccion a la corriente y el pulgar apunta hacia donde va el campo.

Una sola bobina no da un campo homogeneo: el campo es máximo en su centro y se debilita rápido al alejarse.

### 3.2 La configuracion de Helmholtz

Si ponemos dos bobinas iguales, enfrentadas, con la corriente llendo en el mismo sentido, en el espacio entre ellas el campo de una se debilita justo mientras el de la otra se fortalece. A una separación precisa las dos tendencias se compensan y queda una zona central casi plana.

- Para bobinas **circulares**, la separación justa es igual al radio.
- Para bobinas **cuadradas** de lado ``L``, es:

```math
d = 0{,}5445\,L
```

En el grafico se puede obseervaar la relación ``d/L`` donde la franja verde marca dónde el campo se aparta menos de 1 % del valor del centro.
"""

# ╔═╡ 25169c00-231a-41d0-b23c-bfd9db80599e
md"""
Relación de separación ``d/L``: $(@bind r_demo Slider(0.2:0.0005:1.0, default = 0.5445, show_value = true))
"""

# ╔═╡ 660505fb-d06a-44e0-9d96-b5b4993510bb
md"""

- Si las bobinas están **muy juntas**, la suma forma un pico en el centro.
- Si están **muy separadas**, se hunde en el centro y aparecen dos picos.
- En **0,5445** la curva queda lo más plana posible.
"""

# ╔═╡ bd085751-070e-4418-9499-1b19290bfcce
# Relación d/L que usa para la Helmholtz cuadrada
r_helm = 0.5445

# ╔═╡ df2fb722-93d5-4806-bd6a-e841537923c7
md"""
### 3.3 Geometria final de la bobina

Una bobina cuadrada se arma con **cuatro tramos rectos de perfil de aluminio**, que es mucho más fácil de construir que un aro perfecto. Además, tres pares cuadrados se **anidan** fácilmente uno dentro de otro para formar la jaula de 3 ejes. 

### 3.4 Dimencionamiento de las bobinas 

Dentro de la jaula va el giroscopio: la plataforma de pruebas que sostiene al satélite y le permite girar. Este tiene que tener lugar entre las dos bobinas del par más pequeño, para este caso se toma su lado mas extenso de 45cm.

Mirando de lado un par de bobinas: cada bobina está montada en un **perfil** (de 20 mm) centrado en su plano, y queremos una **holgura** (1 cm) entre el giroscopio y el perfil. Desde el centro hasta la bobina (``d/2``) tiene que caber:

```math
\frac{d}{2} \;\ge\; \underbrace{\frac{\text{giroscopio}}{2}}_{22.5} + \underbrace{\text{holgura}}_{1} + \underbrace{\frac{\text{perfil}}{2}}_{1} = 27\ \text{cm} \quad\Longrightarrow\quad d \ge 54\ \text{cm}
```

Y como ``d = 0{,}5445\,L``, el lado mínimo es ``L = d / 0{,}5445``.
"""

# ╔═╡ d110feb3-2397-476c-b92f-0f2925ddacee
md"""
Tamaño del giroscopio (cm): $(@bind g_cm Slider(20:0.5:80, default = 45, show_value = true))

Holgura (cm): $(@bind holgura_cm Slider(0:0.5:5, default = 1, show_value = true))

Ancho/alto del perfil (cm): $(@bind perfil_cm Slider(1:0.5:8, default = 2, show_value = true))
"""

# ╔═╡ e9a3333a-256c-41ef-a501-c960c758e9f3
md""" Asi la separacion de cada par debe ser igual a :"""

# ╔═╡ edfbbd1a-c924-43a1-8088-5d686074f138
begin
	d_min = g_cm + 2holgura_cm + perfil_cm   # cm: separación mínima del par más pequeño
	L_min = d_min / r_helm                   # cm: lado mínimo de ese par
	L_X = ceil(L_min)                        # cm: redondeamos hacia arriba al centímetro
	Δ = 2(perfil_cm + holgura_cm)            # cm: cuánto crece cada par para envolver al anterior
	L_ejes = [L_X, L_X + Δ, L_X + 2Δ]        # cm: lados de los pares X, Y, Z
	d_ejes = r_helm .* L_ejes   # cm: separaciones de los pares X, Y, Z
end

# ╔═╡ 0df888b1-ddc3-4436-b3b0-b20e88063bc6
md"""
Así queda la jaula completa. Cada color es un par de bobinas; la línea punteada que sale del centro indica la dirección del campo que crea ese par.

Mostrar: $(@bind mostrar MultiCheckBox(["Bobinas X", "Bobinas Y", "Bobinas Z", "Giroscopio", "Zona homogénea", "CubeSat"]; default = ["Bobinas X", "Bobinas Y", "Bobinas Z", "Zona homogénea", "CubeSat"]))
"""

# ╔═╡ 4628876c-741e-44c1-99c0-60f2785b6f5a
md"""
---
## 4. ¿Cuánto campo provee una bobina? La ley de Biot–Savart

### 4.1 Principio fisico

La **ley de Biot–Savart** dice cuánto campo magnético produce un diferencial de alambre con corriente:

```math
d\vec B = \frac{\mu_0}{4\pi}\,\frac{I\; d\vec{l} \times \hat r}{r^2}
```


- ``d\vec l``: diferencial de alambre, apuntando en el sentido de la corriente.
- ``r``: la distancia desde el diferencial hasta el punto donde medimos; ``\hat r`` es la dirección.
- El campo es **proporcional a la corriente** ``I``.
- Inverso al cuadrado de la distancia (``1/r^2``}.
- El producto cruz ``\times`` significa que el campo es **perpendicular** tanto al alambre como a la línea que lo une con el punto: da vueltas alrededor del alambre (regla de la mano derecha).

Para obtener el campo de una bobina entera, se suman todas las contribuciones.

### 4.2 Una espira cuadrada, sobre su eje

Si la espira cuadrada tiene lado ``L = 2a`` y medimos sobre el eje que pasa por su centro, a una distancia ``z``, la suma da (para ``N`` vueltas):

```math
B_\text{espira}(z) = \frac{2\mu_0 N I\, a^2}{\pi\,\left(a^2 + z^2\right)\sqrt{2a^2 + z^2}}
```
"""

# ╔═╡ 1fda2d6f-75f7-475c-b099-c13766881b98
md"""
Campo (en T) sobre el eje de UNA espira cuadrada de lado `L` (m), a una distancia `z` (m)
de su centro, por cada amperio-vuelta (NI = 1).
"""

# ╔═╡ fcb17b35-d8fd-4713-acb9-a5c1fe5d67dd
μ0 = 4π * 1e-7   # T·m/A, permeabilidad del vacío

# ╔═╡ 29c6db60-1ed6-4c2c-b2b4-32905282cd3e
B_espira(z, L) = (a = L / 2; 2μ0 * a^2 / (π * (a^2 + z^2) * sqrt(2a^2 + z^2)))

# ╔═╡ 7fbe4978-87c0-453b-a8e5-61f67ffa9597
r_optimo = let
	# Buscamos d/L tal que la 2.ª derivada del campo de UNA espira (L = 1 m)
	# sea cero en z = d/2. Derivada numérica + método de bisección.
	f(z) = B_espira(z, 1.0)
	f2(z; h = 1e-4) = (f(z + h) - 2f(z) + f(z - h)) / h^2
	lo, hi = 0.3, 0.8                 # intervalo donde sabemos que está la respuesta
	for _ in 1:60
		m = (lo + hi) / 2
		f2(lo / 2) * f2(m / 2) <= 0 ? (hi = m) : (lo = m)
	end
	(lo + hi) / 2
end

# ╔═╡ 51a44e7e-2a30-4a1f-b818-912bb15f218c
"""
Campo (en T) sobre el eje de un PAR de espiras cuadradas de lado `L` (m), separadas
`d` (m), con `NI` amperio-vueltas. `z` se mide desde el punto medio entre las dos.
"""
B_par(z, L, d, NI) = NI * (B_espira(z - d / 2, L) + B_espira(z + d / 2, L))

# ╔═╡ f89b9001-46c7-4fdb-b810-117b342059b6
md"""
### 4.3 El par de Helmholtz

Con dos bobinas en ``z = +h/2`` y ``z = -h/2`` (``h`` es la separación, la misma ``d`` de la sección 3) se suman los dos campos. 

```math
\begin{aligned}
B(z) = \frac{2\mu_0 N I a^2}{\pi}\Bigg[ &\frac{1}{\left(a^2+(z-\tfrac{h}{2})^2\right)\sqrt{2a^2+(z-\tfrac{h}{2})^2}} \\[4pt]
+ &\frac{1}{\left(a^2+(z+\tfrac{h}{2})^2\right)\sqrt{2a^2+(z+\tfrac{h}{2})^2}}\Bigg], \qquad L = 2a
\end{aligned}
```

La primera fracción es la bobina de la derecha y la segunda la de la izquierda: cada una es la fórmula de una sola espira, evaluada a la distancia correspondiente.

### 4.4 En el centro (``z = 0``)

En el centro las dos bobinas aportan lo mismo, y todo lo que depende de la geometría se puede agrupar en un solo número ``C``:

```math
B(0) = \mu_0\,N I\cdot C, \qquad C = \frac{4a^2}{\pi\,\left(a^2+\tfrac{h^2}{4}\right)\sqrt{2a^2+\tfrac{h^2}{4}}} \quad [\text{1/m}]
```

``C`` solo depende del tamaño de la bobina.
"""

# ╔═╡ 0de38e59-3c0f-47e8-b716-09749ecee652
# Coeficiente geométrico C (1/m) de un par de lado L (m) con separación d = r·L
C_coef(L; r = r_helm) = B_par(0.0, L, r * L, 1.0) / μ0

# ╔═╡ 5edc19be-f1cd-4de1-a6f8-204f033dba73
md"""

Se observa que ``C`` baja cuando la bobina crece, es decir:

```math
B(0) \propto \frac{NI}{L}
```

Una bobina más grande tiene el alambre más lejos del centro, así que con la misma corriente da **menos campo**. Por eso los pares Y y Z, que son más grandes, necesitarán más corriente que el X.
"""

# ╔═╡ cd93dea0-a34b-4e07-ba27-0a0a13b5ac8c
md"""
---
## 5. Devanado: vueltas(N) y corriente(I)

### 5.1 Amperios-vuelta

De ``B(0) = \mu_0 NI\cdot C`` despejamos lo que necesitamos:

```math
NI = \frac{B_\text{obj}}{\mu_0\, C}
```

``NI`` se mide en **amperios-vuelta (A·vuelta)**. Lo que importa para el campo es el **producto** ``N\times I`` 45 vueltas con 1,64 A dan el mismo campo que 1 vuelta con 73,7 A, o que 10 vueltas con 7,37 A.

Con el campo objetivo ``B_\text{obj}`` de la sección 2 y los lados de la sección 3.Se hace el calculo de los amperios-vuelta que necesita cada eje:
"""

# ╔═╡ 1660aadd-b037-4477-b3cd-0b28c616d0f8
# Amperios-vuelta necesarios en cada par (X, Y, Z) para lograr B_obj en el centro
NI_ejes = [B_obj * 1e-6 / (μ0 * C_coef(L / 100)) for L in L_ejes]

# ╔═╡ 5c827fe0-e03a-4440-9107-ef6b2d0fcab4
md"""
### 5.2 Escoger N y obtener la corriente

**¿Cómo se escoge N?** 
- **Más vueltas** → menos corriente (fuentes más pequeñas y baratas), pero más alambre: más peso, más resistencia eléctrica y más trabajo de bobinado.
- **Menos vueltas** → poco alambre, pero corrientes grandes que calientan los cables y exigen fuentes potentes.

Tambien se toma en cuenta el hecho de que las dos bobinas de cada par se conectan en serie.Se considera que una fuente de corriente maxima de 2A es lo ideal, por tanto se fija N = 50 vueltas por bobina.

Número de vueltas por bobina, ``N``: $(@bind N_vueltas Slider(5:1:200, default = 50, show_value = true))

Asi la corriente necesaria en cada par es:
"""

# ╔═╡ e2b5b115-6061-49ec-ae58-32c77c6ccb10
# Corriente necesaria en cada par (A)
I_ejes = NI_ejes ./ N_vueltas

# ╔═╡ 9dc5475b-fb49-478d-9477-b019e9eb0229
md"""
### 5.3 Alambre necesario


Cada vuelta mide ``4L`` (el perímetro del cuadrado), y cada par tiene ``2N`` vueltas. Con la resistividad del cobre (``\rho = 1{,}68\times10^{-8}\ \Omega\cdot\text{m}``) se obtiene la resistencia ``R``, de ahí el voltaje que debe dar la fuente (``V = IR``) y la potencia que se convierte en calor (``P = I^2R``).

Calibre del alambre: $(@bind awg Select([12 => "AWG 12 (2,05 mm)", 14 => "AWG 14 (1,63 mm)", 16 => "AWG 16 (1,29 mm)", 18 => "AWG 18 (1,02 mm)", 20 => "AWG 20 (0,81 mm)", 22 => "AWG 22 (0,64 mm)"]; default = 20))

asi la masa de alambre para cada eje queda como:
"""

# ╔═╡ 1371faf2-3a2d-459b-8357-af1076cae8c3
begin
	ρ_cu = 1.68e-8                                  # Ω·m, resistividad del cobre a 20 °C
	dens_cu = 8960                                  # kg/m³
	diam_awg(n) = 0.127e-3 * 92^((36 - n) / 39)      # m, fórmula estándar del calibre AWG
	A_alambre = π * diam_awg(awg)^2 / 4              # m², sección del alambre

	largo_ejes = [2 * N_vueltas * 4 * L / 100 for L in L_ejes]   # m de alambre por par
	R_ejes = ρ_cu .* largo_ejes ./ A_alambre                     # Ω
	V_ejes = I_ejes .* R_ejes                                     # V
	P_ejes = I_ejes .^ 2 .* R_ejes                                # W
	masa_ejes = dens_cu .* A_alambre .* largo_ejes               # kg
end

# ╔═╡ 497ac44c-20ce-41ac-9e1a-930c011639fb
md"""
### 5.4 Restriccion de espacio 

El perfil de aluminio de 20 × 20 mm es de **ranura en T**: cada cara tiene un carril con una **boca** angosta entre dos **labios**, y debajo de los labios una **cámara** más ancha. El alambre se enrolla **dentro del carril** de la cara exterior de cada marco: así queda protegido y no se corre de su sitio.

Medidas del carril (en mm):

Ancho de la boca (abertura exterior): $(@bind boca_mm NumberField(1:0.1:20, default = 6))

Ancho de la cámara interior: $(@bind camara_mm NumberField(1:0.1:20, default = 11))

Profundidad total de la ranura, medida desde la cara exterior: $(@bind prof_mm NumberField(1:0.1:20, default =5))

Espesor del labio exterior: $(@bind labio_mm NumberField(0.5:0.1:10, default = 1.5))

De esta forma las vueltas que le caben al perfil de aluminio son :

"""

# ╔═╡ 7b111d0e-e224-44e7-a660-3689a360ba09
begin
	# Carril en T de la cara exterior del perfil (m)
	b_boca = boca_mm / 1000              # ancho de la boca
	b_camara = camara_mm / 1000          # ancho de la cámara
	prof_ranura = prof_mm / 1000         # profundidad total, desde la cara exterior
	e_labio = labio_mm / 1000            # espesor del labio
	h_camara = prof_ranura - e_labio     # alto de la cámara, bajo los labios
	d_esmalte = 1.06 * diam_awg(awg)     # m, diámetro con esmalte grado 2
	ancho_perfil = perfil_cm / 100       # m

	carril = (bcam = b_camara, hcam = h_camara, bboca = b_boca, prof = prof_ranura, ancho = ancho_perfil)
	cara = (bcam = 0.0, hcam = 0.0, bboca = 0.0, prof = 0.0, ancho = ancho_perfil)   # sin carril, para comparar

	"""
	Capas `(altura de su base, espiras)` de `N` espiras de diámetro `d`, llenando desde el fondo
	de la ranura `r`. Las alturas se miden desde el fondo; la cara del perfil está en `r.prof`.
	"""
	function acomodar(N, d; bcam, hcam, bboca, prof, ancho)
		capas = Tuple{Float64,Int}[]
		y = 0.0
		while N > 0
			# el ancho más angosto que toca la capa entre y y y + d
			libre = y + d <= hcam + 1e-9 ? bcam : y >= prof - 1e-9 ? ancho : bboca
			n = min(N, floor(Int, libre / d + 1e-9))
			n == 0 && error("El alambre de $(round(1e3 * d, digits = 2)) mm no pasa por la boca del carril")
			push!(capas, (y, n))
			N -= n
			y += d
		end
		capas
	end

	"Centros `(x, y)` (m) de las espiras de cada capa, centradas en el carril."
	centros(capas, d) = [((i - (n + 1) / 2) * d, y + d / 2) for (y, n) in capas for i in 1:n]

	"Espiras de diámetro `d` que caben enteras dentro de la ranura `r`."
	vueltas_en(r, d) = sum((n for (y, n) in acomodar(10^4, d; r...) if y + d <= r.prof + 1e-9); init = 0)

	"Distancia media geométrica (m) de un haz de alambres de radio de cobre `rcu` con centros `c`."
	function dmg(c, rcu)
		s = length(c) * log(exp(-1 / 4) * rcu)   # cada alambre consigo mismo
		for i in eachindex(c), j in eachindex(c)
			i != j && (s += log(hypot(c[i][1] - c[j][1], c[i][2] - c[j][2])))
		end
		exp(s / length(c)^2)
	end

	capas_carril = acomodar(N_vueltas, d_esmalte; carril...)
	capas_cara = acomodar(N_vueltas, d_esmalte; cara...)
	N_carril = vueltas_en(carril, d_esmalte)                     # máximo que cabe en el carril
	sobresale = max(0.0, capas_carril[end][1] + d_esmalte - prof_ranura)   # m, alto del haz por fuera de la cara
	DMG = dmg(centros(capas_carril, d_esmalte), diam_awg(awg) / 2)          # m
	DMG_cara = dmg(centros(capas_cara, d_esmalte), diam_awg(awg) / 2)       # m, si fuera sobre la cara
	N_carril
end

# ╔═╡ 90093a54-0f76-4ea2-8ad3-457fc5f3c93f
md"""
### 5.5 Inductancia de cada par

La inductancia mide cuánto se opone una bobina a que su corriente cambie. Cuando la corriente sube o baja, el campo que guarda la bobina también cambia, y ese cambio induce un voltaje que va en contra. Importa en dos cosas:
- la corriente no llega a su valor de inmediato: tarda un tiempo del orden de ``\tau = L_\text{par}/R``, la constante de tiempo del par (``R`` es la resistencia de la sección 5.3);
- mientras la corriente cambia, la fuente debe dar un voltaje extra ``L_\text{par}\,dI/dt``.

Las dos bobinas de un par están en serie y con la corriente en el mismo sentido, así que el campo de cada una también atraviesa a la otra. A la inductancia propia de cada bobina se le suma la inductancia mutua``M`` entre ellas, dos veces:

```math
L_\text{par} = 2\,\left(L_\text{bobina} + M\right)
```

**Inductancia propia.** Para una bobina cuadrada de ``N`` vueltas y lado ``L`` (aquí ``L`` sin subíndice sigue siendo el lado):

```math
L_\text{bobina} = N^2\,\frac{2\mu_0 L}{\pi}\left[\ln\frac{L}{\text{DMG}} - 0{,}774\right]
```

La fórmula sale de sumar las inductancias de los 4 lados rectos y restar la mutua entre lados opuestos, por los que la corriente va en sentidos contrarios.

DMG es la **distancia media geométrica** de la sección del haz de alambre: un promedio de las distancias entre todos los pares de alambres del haz. Se calcula con las posiciones de las espiras en el carril (sección 5.4):

```math
\ln \text{DMG} = \frac{1}{N^2}\sum_{i=1}^{N}\sum_{j=1}^{N} \ln r_{ij}
```

``r_{ij}`` es la distancia entre los centros de los alambres ``i`` y ``j``. Cada alambre consigo mismo cuenta como ``r_{ii} = 0{,}7788`` veces su radio de cobre. Un haz más compacto tiene una DMG más pequeña y por eso un poco más de inductancia.

**Inductancia mutua.** Cada lado de una bobina tiene en la otra un lado paralelo a distancia ``d`` con la corriente en el mismo sentido (suma) y otro a distancia ``\sqrt{L^2 + d^2}`` con la corriente al revés (resta). Los lados perpendiculares no se acoplan. Con la mutua de dos tramos rectos paralelos de largo ``l`` separados una distancia ``r``,

```math
m(l, r) = \frac{\mu_0}{2\pi}\left[\,l\,\ln\frac{l + \sqrt{l^2 + r^2}}{r} - \sqrt{l^2 + r^2} + r\right]
```

queda

```math
M = 4N^2\left[\,m(L, d) - m\!\left(L, \sqrt{L^2 + d^2}\right)\right]
```
"""

# ╔═╡ 32f8e23c-40ed-4653-825d-66d13c109f48
begin
	# Inductancia mutua (H) entre dos tramos rectos paralelos de largo l, enfrentados a distancia r
	m_tramos(l, r) = μ0 / 2π * (l * log((l + hypot(l, r)) / r) - hypot(l, r) + r)

	Lbob_ejes = [N_vueltas^2 * 2μ0 * (L / 100) / π * (log(L / 100 / DMG) - 0.774) for L in L_ejes]   # H, una bobina
	M_ejes = [4N_vueltas^2 * (m_tramos(L / 100, d / 100) - m_tramos(L / 100, hypot(L, d) / 100))
	          for (L, d) in zip(L_ejes, d_ejes)]                                                    # H, entre las dos bobinas
	Lpar_ejes = 2 .* (Lbob_ejes .+ M_ejes)   # H, par en serie
	τ_ejes = Lpar_ejes ./ R_ejes             # s, constante de tiempo
end

# ╔═╡ 478ea9ab-a3f2-472a-8d97-e4d3684d218b
md"""
### 5.6 Calentamiento: disipacion  del devanado

Toda la potencia eléctrica que entra a las bobinas se convierte en **calor** (efecto Joule). Un campo magnético que no cambia no consume energía, así que nada de esa potencia se aprovecha. El calor sube la temperatura del alambre y del perfil hasta llegar a un equilibrio: la temperatura en la que la bobina le entrega al aire exactamente el mismo calor que genera.

**Calor generado.** Es la potencia ``P = I^2 R`` de la sección 5.3, agregando que la resistencia del cobre crece con la temperatura,

```math
R(T) = R_{20}\left[1 + \alpha\,\left(T - 20\,^\circ\text{C}\right)\right], \qquad \alpha_\text{Cu} = 0{,}00404\ \text{K}^{-1}
```

``R_{20}`` es la resistencia a 20 °C (la de la tabla de la sección 5.3). La corriente la fija la fuente, así que la bobina caliente genera un poco más de calor que la fría.

**Calor que sale.** La bobina pierde calor por su superficie exterior ``A`` de dos formas:
- **convección natural**: el aire que toca la bobina se calienta, sube y deja entrar aire más frío;
- **radiación**: toda superficie caliente emite radiación infrarroja.

Las dos crecen con la diferencia de temperatura entre la bobina y el aire del laboratorio:

```math
Q = \left(h_\text{conv} + h_\text{rad}\right)A\,\left(T - T_\text{amb}\right)
```

``h`` es el **coeficiente de transferencia de calor**, en W/(m²·K): cuántos vatios salen por cada metro cuadrado de superficie y por cada grado de diferencia. La temperatura de equilibrio es la ``T`` en la que ``P(T) = Q(T)``, y se busca con el mismo método de bisección de la sección 3.

**Superficie.** El haz de alambre va dentro del carril (sección 5.4) y le pasa su calor al aluminio por el fondo, los lados y los labios. El aluminio conduce muy bien el calor, así que toda la sección (perfil + haz) se calienta casi por igual y entrega calor al aire por su borde exterior. Ese borde es el contorno del perfil (4 veces su ancho) más los dos lados de la parte del haz que sobresalga de la cara. La superficie de una bobina es ``A = 4L\times`` (borde de la sección).
"""

# ╔═╡ 4179b755-8888-4ac1-9e90-85e98b1f90f4
md"""
Temperatura ambiente del laboratorio (°C): $(@bind T_amb Slider(10:1:40, default = 20, show_value = true))

Acabado del perfil: $(@bind ε_perfil Select([0.8 => "Anodizado (ε ≈ 0,8)", 0.1 => "Aluminio sin anodizar (ε ≈ 0,1)", 0.95 => "Pintado de negro mate (ε ≈ 0,95)"]; default = 0.8))
"""

# ╔═╡ cb618cf5-b404-4d3f-af27-b74bfc1fcee8
md"""

---
## 6. Comprobación: Uniformidad del campo

La fórmula de la sección 4 solo da el campo sobre el eje (una línea). Pero el satélite ocupa un **volumen**, así que hay que mirar también los puntos fuera del eje.

Para eso usamos otra vez Biot–Savart, ahora para un **tramo recto** visto desde **cualquier punto** del espacio.Cada bobina cuadrada son 4 tramos, y sumando los 8 tramos de un par obtenemos el campo en cualquier punto ``(x, y, z)``.
"""

# ╔═╡ 21cb32e6-3fb5-4475-8e20-575796d24dfe
"Campo (T) en el punto `P` debido a un tramo recto de `A` a `B` con corriente `I`."
function B_tramo(P, A, B, I)
	a, b = A - P, B - P
	na, nb = norm(a), norm(b)
	μ0 * I / 4π * (na + nb) / (na * nb * (na * nb + dot(a, b))) * cross(a, b)
end

# ╔═╡ 24b6b149-97cf-4d2a-8605-5fbe98a9cf29
"""
Campo (T) en `P` de un par de espiras cuadradas de lado `L` (m), separación `d` (m)
y `NI` amperio-vueltas. El eje del par es el eje z del sistema local.
"""
function B_par3D(P, L, d, NI)
	h = L / 2
	B = zeros(3)
	for z0 in (-d / 2, d / 2)
		esquinas = ([h, h, z0], [-h, h, z0], [-h, -h, z0], [h, -h, z0])
		for i in 1:4
			B += B_tramo(P, esquinas[i], esquinas[mod1(i + 1, 4)], NI)
		end
	end
	B
end

# ╔═╡ 76baad62-fa2f-4540-9129-bded98520c8d
md"""
Escoge qué par revisar y con qué separación. El mapa es un corte que contiene el eje del par y muestra cuánto se aparta el campo del valor que tiene en el centro:
- **blanco**: igual que en el centro,
- **rojo**: más fuerte, **azul**: más débil (la escala se satura en ±2 %),
- las **líneas negras** marcan dónde la diferencia llega a 1 %,
- el **círculo verde** es la zona homogénea de la sección 1, donde va el satélite.
"""

# ╔═╡ 75724c50-ee26-4c46-98e9-8c12ed7110fd
md"""
Par: $(@bind eje_ver Select(["X", "Y", "Z"]))

Relación ``d/L``: $(@bind r_ver Slider(0.35:0.0005:0.8, default = 0.5445, show_value = true))
"""

# ╔═╡ 025766ec-d4c4-4673-a31c-468771a2fec0
begin
	k_ver = findfirst(==(eje_ver), ["X", "Y", "Z"])
	L_ver = L_ejes[k_ver] / 100              # m
	d_ver = r_ver * L_ver                    # m
	B_centro = B_par3D([0.0, 0, 0], L_ver, d_ver, 1.0)

	# Desviación de magnitud (%, con signo) y de dirección (°) respecto al centro
	function desviacion(P)
		B = B_par3D(P, L_ver, d_ver, 1.0)
		100 * (norm(B) / norm(B_centro) - 1),
		acosd(clamp(dot(B, B_centro) / (norm(B) * norm(B_centro)), -1, 1))
	end

	# Peor caso dentro de la esfera de diámetro D (puntos de una rejilla 3D)
	R_esf = D_cm / 200
	pts = [[x, y, z] for x in range(-R_esf, R_esf, length = 17),
	                     y in range(-R_esf, R_esf, length = 17),
	                     z in range(-R_esf, R_esf, length = 17) if x^2 + y^2 + z^2 <= R_esf^2]
	devs = desviacion.(pts)
	dev_max_mag = maximum(abs ∘ first, devs) # %
	dev_max_ang = maximum(last, devs)        # grados
end

# ╔═╡ 0dda238e-8351-4e42-a4c9-c18c5a36a59c
md"""
**¿Por qué el mapa tiene esa forma de flor?** Con la separación de Helmholtz, al alejarse del centro el campo se debilita si uno se mueve a lo largo del eje o hacia los lados (zonas azules), y se fortalece  si se mueve en diagonal (zonas rojas). En la frontera entre unas y otras el campo coincide con el del centro, y por eso el blanco se estira en forma de X. Lo importante es que cerca del centro todo es casi blanco: ahí va el satélite.

Un criterio habitual en jaulas de Helmholtz es que el campo varíe menos de ~1 % en el volumen de prueba. Con la geometría de el CubeSat 1U queda **muy por debajo** de ese límite, porque la jaula (≈ 1 m) es mucho más grande que el satélite (≈ 18 cm).
"""

# ╔═╡ 08c18c59-d31b-4aa5-b61b-d9206b1fb174
md"""
---
## 7. Restriccion: nada ferromagnético

la estructura no puede tener materiales ferromagnéticos.

Los materiales ferromagnéticos (hierro, acero al carbono, níquel, cobalto) se imantan cuando están dentro de un campo magnético, y entonces producen su propio campo. Una pieza de hierro cerca del centro concentra y tuerce las líneas de campo, y arruina justamente la uniformidad que tanto cuidamos.

|Se puede usar |  Evitar |
|:---|:---|
| Aluminio (perfiles) | Hierro, acero al carbono |
| Plásticos, impresiones 3D, madera | Tornillos y tuercas de acero comunes |
| Latón, bronce, tornillos de nylon | Piezas de níquel o niqueladas |
| Cobre (el alambre) | Herramientas y muebles metálicos cerca de la jaula |

Tambien algunos inoxidables (los ferríticos, y a veces el 304 después de doblado o mecanizado) son algo magnéticos.
"""

# ╔═╡ 69d3dc2e-1904-44dc-9207-a7ae8c7e995d
md"""
---
## 8. Resumen del diseño

Esta tabla toma los valores de **todos los controles de arriba**: si cambiaste algo, aquí aparece reflejado.
"""

# ╔═╡ 2fa60789-9ea4-471b-99aa-2b69584bb79e
md"""
---
## Apéndice: herramientas del cuaderno

Estas celdas definen funciones auxiliares (formato de números, figuras 3D, desplegables y gráficas que funcionan sin internet). 
"""

# ╔═╡ 127d07c1-71d9-478f-b254-6398e757c941
html"""
<style>
/* Tablas con la misma letra que el texto (Pluto usa letra de código por defecto) */
pluto-output .markdown table > tbody td { font-family: inherit; font-size: 0.95em; font-variant-ligatures: normal; }
</style>
"""

# ╔═╡ ea853ecb-7073-4ffe-9b03-6b4d40d7b0ec
"Formatea un número con `n` decimales."
fmt(x, n = 2) = Printf.format(Printf.Format("%.$(n)f"), x)

# ╔═╡ 734d8df8-56ca-4d4c-a3ef-09e737943bde
Markdown.parse("""
!!! tip "Resultado"
    Diámetro de la zona homogénea necesaria: **D = $(fmt(D_cm)) cm**

En la figura se observa la caja del satélite, su diagonal en verde y la esfera que la envuelve.
""")

# ╔═╡ 7de83929-3370-4c88-99ab-06dfc6f9ebb6
Markdown.parse("""
Con el ángulo elegido (``\\theta`` = $(θ_deg)°), las bobinas deben generar **$(fmt(B_bobinas_θ, 1)) µT**.

El **peor caso** es ``\\theta = 180°``: el campo deseado apunta justo al revés del de Bogotá. Las bobinas tienen que anular los $(B_tierra) µT de la Tierra **y** además crear los $(B_orbita) µT de la órbita:

```math
B_\\text{peor} = B_\\text{Tierra} + B_\\text{órbita} = $(B_tierra) + $(B_orbita) = $(fmt(B_peor, 0))\\ \\mu\\text{T}
```

Se diseña con margen, en este caso se elige un factor de seguridad de 2, puesto que el campo real del laboratorio varía (equipos, estructuras metálicas, cables cercanos), los cálculos son aproximados, conviene poder simular campos más fuertes que el nominal y no es bueno operar las bobinas al límite dada la disipacion de calor.

!!! tip "Campo objetivo de diseño"
    ``B_\\text{obj}`` = $(fmt(FS, 1)) × $(fmt(B_peor, 0)) = **$(fmt(B_obj, 0)) µT por eje**
""")

# ╔═╡ 517de74b-0cf1-4f0b-8e82-643c2af149e5
Markdown.parse("""
**Par más pequeño (eje X):**

```math
d_\\text{mín} = $(fmt(g_cm, 1)) + 2\\cdot $(fmt(holgura_cm, 1)) + $(fmt(perfil_cm, 1)) = $(fmt(d_min, 1))\\ \\text{cm}
\\qquad
L_\\text{mín} = \\frac{$(fmt(d_min, 1))}{0.5445} = $(fmt(L_min, 2))\\ \\text{cm}
```

!!! tip "Lado del par X"
    ``L_\\text{mín}`` = $(fmt(L_min, 2)) cm, que se redondea **hacia arriba** a **$(fmt(L_X, 0)) cm**.

**Los otros dos pares.** Cada par tiene que rodear el anterior. En cada lado debe dejar espacio para el perfil del par interior más la holgura, así que el lado crece

```math
\\Delta = 2\\,(\\text{perfil} + \\text{holgura}) = 2\\,($(fmt(perfil_cm, 1)) + $(fmt(holgura_cm, 1))) = $(fmt(Δ, 1))\\ \\text{cm}
```
""")

# ╔═╡ 5da0dbdb-5db2-4b1c-bfa4-f841ea69f5eb
let
	todas = [n for (y, n) in acomodar(10^4, d_esmalte; carril...) if y + d_esmalte <= prof_ranura + 1e-9]
	estado = if N_vueltas == N_carril
		"Las **$(N_vueltas) vueltas caben** justas dentro del carril: lo llenan por completo."
	elseif N_vueltas < N_carril
		"Las **$(N_vueltas) vueltas caben** completas dentro del carril, y todavía quedaría lugar para $(N_carril - N_vueltas) más."
	else
		"**No caben las $(N_vueltas) vueltas.** Quedan $(N_carril) dentro del carril y las otras $(N_vueltas - N_carril) asoman por la boca o van encima de la cara, así que el haz sobresale $(fmt(1e3 * sobresale, 1)) mm de la cara del perfil. Para que quepan todas hay que bajar a $(N_carril) vueltas, con más corriente, o usar un alambre más delgado."
	end
	Markdown.parse("""
	!!! tip "Devanado de AWG $(awg) en el carril"
	    En el carril caben **$(N_carril) vueltas** de AWG $(awg) (alambre de $(fmt(1e3 * d_esmalte, 2)) mm con esmalte), en capas de $(join(todas, ", ")) espiras de abajo hacia arriba.

	$(estado)

	""")
end

# ╔═╡ 51ec3d44-bc87-4783-85f5-afa0190995b0
Markdown.parse("""
!!! tip "Devanado de $(N_vueltas) vueltas de AWG $(awg)"
    Con las espiras en el carril, el haz tiene una DMG de **$(fmt(1e3 * DMG, 2)) mm**. Sobre la cara plana del perfil sería de $(fmt(1e3 * DMG_cara, 2)) mm: en el carril el haz es más compacto y la inductancia propia de cada bobina sale un $(fmt(100 * (log(L_ejes[1] / 100 / DMG) - 0.774) / (log(L_ejes[1] / 100 / DMG_cara) - 0.774) - 100, 0)) % mayor.

La inductancia crece con ``N^2``: el doble de vueltas da casi cuatro veces más inductancia y el doble de resistencia, así que ``\\tau`` también se duplica. La mutua aporta poco, alrededor del $(fmt(100 * maximum(M_ejes ./ (Lbob_ejes .+ M_ejes)), 0)) % del total, porque las dos bobinas están separadas más de medio lado.

Con ``\\tau`` de a lo sumo $(fmt(1e3 * maximum(τ_ejes), 1)) ms, la corriente llega al 99 % de su valor en ``5\\tau \\approx`` $(fmt(5e3 * maximum(τ_ejes), 0)) ms. Es mucho más rápido que los cambios de campo que ve un CubeSat al girar en órbita, así que la inductancia no limita las pruebas.
""")

# ╔═╡ 68f7d89d-b0ef-4a5c-88f2-e5b17a7d96b3
Markdown.parse("""
Las barras rojas son las dos bobinas vistas de canto.

!!! tip "Dentro de la esfera de $(fmt(D_cm, 1)) cm, par $eje_ver"
    - el **tamaño** del campo varía como máximo **$(fmt(dev_max_mag, 3)) %**
    - la **dirección** del campo se desvía como máximo **$(fmt(dev_max_ang, 3))°**
""")

# ╔═╡ 96a3a5aa-183c-41c4-87d8-dc75d029e895
"Construye una tabla Markdown a partir de un encabezado y una lista de filas."
function tabla(encabezado, filas)
	s = "| " * join(encabezado, " | ") * " |\n|" * join(fill(":---", length(encabezado)), "|") * "|\n"
	for f in filas
		s *= "| " * join(string.(f), " | ") * " |\n"
	end
	Markdown.parse(s)
end

# ╔═╡ b03602d8-6c46-4061-9d46-d4cefa05a519
tabla(["Eje", "Lado ``L``", "Separación ``d = 0{,}5445\\,L``"],
      [[eje, fmt(L / 100, 2) * " m", fmt(d, 2) * " cm"]
       for (eje, L, d) in zip(["X", "Y", "Z"], L_ejes, d_ejes)])

# ╔═╡ 5db17c51-8b34-4670-b558-3da87d4d4e7d
let
	Ls = [0.9, 0.96, 1.02]
	notas = [1.29583, 1.19985, 1.11710]
	tabla(["Lado  `` L ``", "`` C\\cdot\\pi ``", "``  C``"],
	      [[fmt(L, 2) * " m", fmt(C_coef(L) * π, 4), fmt(C_coef(L), 5)] for (L, n) in zip(Ls, notas)])
end

# ╔═╡ a7b47a0a-22b3-447e-885c-4b61696c72ff
tabla(["Eje", "Lado ``L``", "``C`` (1/m)", "``NI`` requerido", "Corriente ``I`` con N = $(N_vueltas)"],
      [[eje, fmt(L / 100, 2) * " m", fmt(C_coef(L / 100), 5), fmt(NI, 2) * " A·vuelta", fmt(I, 2) * " A"]
       for (eje, L, NI, I) in zip(["X", "Y", "Z"], L_ejes, NI_ejes, I_ejes)])

# ╔═╡ 53e661ca-0e66-4f25-a6bb-258ce8a758b9
tabla(["Eje", "Alambre por par", "Resistencia", "Voltaje", "Potencia", "Masa de cobre"],
      [[eje, fmt(l, 0) * " m", fmt(R, 1) * " Ω", fmt(V, 1) * " V", fmt(P, 1) * " W", fmt(m, 2) * " kg"]
       for (eje, l, R, V, P, m) in zip(["X", "Y", "Z"], largo_ejes, R_ejes, V_ejes, P_ejes, masa_ejes)])

# ╔═╡ c7bb36c0-a2ac-42a0-a1ed-a11b6dd8bfdb
tabla(["Eje", "Inductancia de una bobina", "Mutua ``M``", "Inductancia del par ``L_\\text{par}``", "Constante de tiempo ``\\tau = L_\\text{par}/R``"],
      [[eje, fmt(1e3 * Lb, 1) * " mH", fmt(1e3 * M, 2) * " mH", fmt(1e3 * Lp, 1) * " mH", fmt(1e3 * τ, 2) * " ms"]
       for (eje, Lb, M, Lp, τ) in zip(["X", "Y", "Z"], Lbob_ejes, M_ejes, Lpar_ejes, τ_ejes)])

# ╔═╡ 70e73562-6b5f-468d-b495-0ffd677983f1
"Semiancho (m) de la zona del eje donde el campo se aparta menos de `tol` del centro (L = 1 m)."
function semiancho_uniforme(r; L = 1.0, tol = 0.01)
	B0 = B_par(0.0, L, r * L, 1.0)
	for z in range(0, L, length = 4001)
		abs(B_par(z, L, r * L, 1.0) / B0 - 1) > tol && return z
	end
	L
end

# ╔═╡ 45792165-6345-439b-ae15-d07a82184fd4
Markdown.parse("""
Con ``d/L`` = $(r_demo), la zona donde el campo (sobre el eje) se aparta menos de 1 % mide **$(fmt(200semiancho_uniforme(r_demo), 1)) cm** de largo. Las líneas punteadas verticales marcan dónde están las bobinas.
""")

# ╔═╡ 2f65f759-0037-41f9-8eec-c4fb10dda0ed
"""
Resistencia térmica (K por cada W/m de bobina) entre el punto más caliente del haz y el
aluminio, por diferencias finitas en la sección del carril `r` con las espiras de `capas`.
La sección se divide en cuadritos de lado `h`. El aluminio queda a temperatura fija y el
aire de afuera se toma aislante, que es el caso más desfavorable.
"""
function resistencia_haz(capas, d, r; k, k_aire, h = 1e-4)
	nx, ny = round(Int, r.ancho / h), ceil(Int, max(r.prof, capas[end][1] + d) / h)
	xc(i) = (i - 0.5) * h - nx * h / 2   # m, centro del cuadrito; x = 0 en el centro del carril
	yc(j) = (j - 0.5) * h                # m, desde el fondo del carril
	alambre = [any(y0 <= yc(j) < y0 + d && abs(xc(i)) < n * d / 2 for (y0, n) in capas) for i in 1:nx, j in 1:ny]
	hueco(i, j) = yc(j) < r.hcam ? abs(xc(i)) < r.bcam / 2 : yc(j) < r.prof && abs(xc(i)) < r.bboca / 2
	# conductividad de cada cuadrito: Inf = aluminio (temperatura fija), 0 = aire de afuera
	kc = [alambre[i, j] ? k : hueco(i, j) ? k_aire : yc(j) < r.prof ? Inf : 0.0 for i in 1:nx, j in 1:ny]
	vecino(i, j) = j < 1 ? Inf : 1 <= i <= nx && j <= ny ? kc[i, j] : 0.0   # bajo el fondo hay aluminio
	libres = [I for I in CartesianIndices(kc) if 0 < kc[I] < Inf]
	q = 1 / count(alambre)               # W/m que genera cada cuadrito de alambre
	T = zeros(nx, ny)
	for _ in 1:100_000                   # Gauss–Seidel con sobrerrelajación
		cambio = 0.0
		for I in libres
			i, j = Tuple(I)
			entra, g_total = alambre[I] ? q : 0.0, 0.0
			for (a, b) in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1))
				kv = vecino(a, b)
				kv == 0 && continue
				g = isinf(kv) ? 2kc[I] : 2kc[I] * kv / (kc[I] + kv)   # conductancia entre los centros
				entra += isinf(kv) ? 0.0 : g * T[a, b]
				g_total += g
			end
			nuevo = entra / g_total
			cambio = max(cambio, abs(nuevo - T[I]))
			T[I] += 1.9 * (nuevo - T[I])
		end
		cambio < 1e-9 && break
	end
	maximum(T[alambre])
end

# ╔═╡ cf499b0f-54fd-4d42-9a03-a53fa848b539
begin
	α_cu = 0.00404                 # 1/K, cuánto sube la resistencia del cobre por cada grado
	c_cu = 385                     # J/(kg·K), calor específico del cobre (alambre)
	c_alum = 897                   # J/(kg·K), calor específico del aluminio (perfil)
	dens_alum = 2700               # kg/m³, densidad del aluminio (perfil)
	σ_SB = 5.670e-8                # W/(m²·K⁴), constante de Stefan–Boltzmann
	p_lab = 75.0                   # kPa, presión del aire en Bogotá (≈ 2600 m)
	ε_esmalte = 0.9                # emisividad del esmalte del alambre
	k_haz = 0.3                    # W/(m·K), conductividad de un haz de alambre esmaltado sin resina
	k_aire = 0.026                 # W/(m·K), aire quieto en los huecos del carril
	lados_verticales = [2, 2, 0]   # los pares X y Y están en planos verticales; el Z, acostado
	A_perfil = 0.4225 * ancho_perfil^2                # m², aluminio de un perfil ranurado tipo 2020

	# Sección de una bobina vista desde afuera: el perfil más lo que sobresale del haz
	function seccion(capas, d, r)
		fuera = max(0.0, capas[end][1] + d - r.prof)   # m, alto del haz por fuera de la cara
		# ancho de alambre que se ve desde afuera: la boca del carril o la capa más ancha que asoma
		visto = maximum((n * d for (y, n) in capas if y + d > r.prof + 1e-9); init = r.bboca)
		borde = 4r.ancho + 2fuera                      # m
		alambre = visto + 2fuera                       # m, parte del borde que es alambre
		(D = r.ancho + fuera,                          # m, alto de la sección
		 borde = borde,
		 ε = (ε_esmalte * alambre + ε_perfil * (borde - alambre)) / borde)   # emisividad promedio
	end
	sec_carril = seccion(capas_carril, d_esmalte, carril)
	sec_cara = seccion(capas_cara, d_esmalte, cara)

	kelvin(T) = T + 273.15
	# Convección natural (W/m²K) de una bobina de lado Lm (m) con nv lados verticales
	h_conv(ΔT, Lm, nv, s) = sqrt(p_lab / 101.325) *
		((4 - nv) * 1.32 * (ΔT / s.D)^0.25 + nv * 1.42 * (ΔT / Lm)^0.25) / 4
	h_rad(T, s) = s.ε * σ_SB * (kelvin(T)^2 + kelvin(T_amb)^2) * (kelvin(T) + kelvin(T_amb))
	calor_generado(T, I, R20) = I^2 * R20 * (1 + α_cu * (T - 20))     # W, una bobina
	calor_perdido(T, Lm, nv, s) = T <= T_amb ? 0.0 :
		(h_conv(T - T_amb, Lm, nv, s) + h_rad(T, s)) * 4Lm * s.borde * (T - T_amb)   # W, una bobina

	# Temperatura de equilibrio de una bobina (°C), por bisección
	function T_equilibrio(I, R20, Lm, nv, s)
		lo, hi = float(T_amb), T_amb + 500.0
		for _ in 1:60
			m = (lo + hi) / 2
			calor_generado(m, I, R20) > calor_perdido(m, Lm, nv, s) ? (lo = m) : (hi = m)
		end
		(lo + hi) / 2
	end

	R20_bob = R_ejes ./ 2                                                    # Ω, una bobina a 20 °C
	T_ejes(s) = [T_equilibrio(I, R, L / 100, nv, s)
	             for (I, R, L, nv) in zip(I_ejes, R20_bob, L_ejes, lados_verticales)]   # °C
	Tfin_ejes = T_ejes(sec_carril)                                           # °C
	Rcal_ejes = R_ejes .* (1 .+ α_cu .* (Tfin_ejes .- 20))                   # Ω, par en equilibrio
	Vcal_ejes = I_ejes .* Rcal_ejes                                          # V
	Pcal_ejes = I_ejes .^ 2 .* Rcal_ejes                                     # W, par en equilibrio
	# Punto más caliente: el calor del haz lo cruza hasta el aluminio (K por cada W/m de bobina)
	Rhaz_carril = resistencia_haz(capas_carril, d_esmalte, carril; k = k_haz, k_aire)
	Tmax_ejes = Tfin_ejes .+ Rhaz_carril .* [P / 2 / (4L / 100) for (P, L) in zip(Pcal_ejes, L_ejes)]   # °C

	# Lo mismo con el haz sobre la cara del perfil, sin carril, para comparar
	Tcara_ejes = T_ejes(sec_cara)
	Pcara_ejes = I_ejes .^ 2 .* R_ejes .* (1 .+ α_cu .* (Tcara_ejes .- 20))
	Rhaz_cara = resistencia_haz(capas_cara, d_esmalte, cara; k = k_haz, k_aire)
	Tmax_cara = Tcara_ejes .+ Rhaz_cara .* [P / 2 / (4L / 100) for (P, L) in zip(Pcara_ejes, L_ejes)]

	# Calentamiento en el tiempo: C dT/dt = P − Q, en pasos de 1 s durante 8 h
	C_ejes = [c_cu * m / 2 + c_alum * dens_alum * A_perfil * 4L / 100
	          for (m, L) in zip(masa_ejes, L_ejes)]                          # J/K, una bobina (alambre + perfil)
	t_min = 0:480                                                            # min
	curvas_T = map(1:3) do k
		Lm, T = L_ejes[k] / 100, float(T_amb)
		curva = [T]
		for _ in 1:last(t_min)
			for _ in 1:60
				T += (calor_generado(T, I_ejes[k], R20_bob[k]) - calor_perdido(T, Lm, lados_verticales[k], sec_carril)) / C_ejes[k]
			end
			push!(curva, T)
		end
		curva
	end
	# Minutos hasta completar el 95 % del calentamiento
	t95_ejes = [t_min[something(findfirst(≥(T_amb + 0.95(Tf - T_amb)), c), length(t_min))]
	            for (c, Tf) in zip(curvas_T, Tfin_ejes)]

	Pcal_ejes
end

# ╔═╡ 034bf3c5-e6d0-4b9c-81c6-40c1b241e749
tabla(["Eje", "Corriente", "Calor del par a 20 °C", "Calor del par en equilibrio", "Temperatura del devanado", "Punto más caliente", "Tiempo al 95 %"],
      [[eje, fmt(I, 2) * " A", fmt(P20, 1) * " W", fmt(P, 1) * " W", fmt(T, 1) * " °C", fmt(Tm, 1) * " °C", "$(t) min"]
       for (eje, I, P20, P, T, Tm, t) in zip(["X", "Y", "Z"], I_ejes, P_ejes, Pcal_ejes, Tfin_ejes, Tmax_ejes, t95_ejes)])

# ╔═╡ 14bb9727-15ca-4954-9de4-76566726a9af
let
	T_peor = maximum(Tmax_ejes)
	veredicto = if T_peor < 45
		"**Es poco.** El esmalte de los alambres para bobinas aguanta en servicio continuo 130 °C o más (según su clase térmica), y una superficie metálica por debajo de unos 45 °C se puede tocar sin problema."
	elseif T_peor < 130
		"**Se puede operar, con cuidado.** El devanado pasa de unos 45 °C, que ya quema al tocarlo por un rato, pero sigue por debajo de los 130 °C que aguanta el esmalte de los alambres para bobinas."
	else
		"**Es demasiado.** El devanado pasaría de 130 °C, el límite de los esmaltes más comunes. Hay que poner más vueltas o un alambre más grueso."
	end
	Markdown.parse("""
	!!! tip "Calor del devanado con el aire a $(T_amb) °C"
	    - Con los tres pares encendidos a la vez a su corriente máxima, la jaula disipa **$(fmt(sum(Pcal_ejes), 0)) W** en total ($(join(fmt.(Pcal_ejes, 1), " / ")) W en los pares X / Y / Z).
	    - Cada devanado se estabiliza entre **$(fmt(minimum(Tfin_ejes), 1)) y $(fmt(maximum(Tfin_ejes), 1)) °C**, es decir, de $(fmt(minimum(Tfin_ejes) - T_amb, 1)) a $(fmt(maximum(Tfin_ejes) - T_amb, 1)) K por encima del aire. El interior del haz llega a lo sumo a $(fmt(T_peor, 1)) °C.
	    - Tarda unos **$(maximum(t95_ejes)) min** en completar el 95 % de ese calentamiento.

	$(veredicto)

	**La resistencia sube con la temperatura.** En equilibrio la resistencia de los pares es $(join(fmt.(Rcal_ejes, 1), " / ")) Ω (a 20 °C era $(join(fmt.(R_ejes, 1), " / ")) Ω), así que la fuente debe poder dar hasta **$(fmt(maximum(Vcal_ejes), 1)) V**. Conviene usar la fuente en **modo de corriente constante (CC)**: el campo depende solo de la corriente, y así no cambia mientras el devanado se calienta. Con voltaje constante, la corriente y el campo bajarían a medida que sube la resistencia.

	De ser necesario debemos estar al tanto de :¿Cómo se reduce el calor? si reemplazamos ``I = NI/N`` y la resistencia de la sección 5.3 en ``P = I^2R``:

	```math
	P_\\text{par} = \\frac{(NI)^2\\; 8\\,\\rho\\,L}{N\\,A_\\text{alambre}}
	```

	``NI`` lo fija el campo (sección 5.1), así que el calor solo depende de ``N\\,A_\\text{alambre}``: el área total de cobre en la sección del devanado. Duplicar las vueltas, o usar un alambre del doble de sección, reduce el calor a la mitad, a cambio del doble de masa de alambre.
	""")
end

# ╔═╡ 845a1707-f08a-43b5-a5e8-ac71114ec5cd
tabla(["Parámetro", "Valor", "De dónde sale"],
      [["Medidas del CubeSat", "$(fmt(a_cm, 2)) × $(fmt(b_cm, 2)) × $(fmt(c_cm, 2)) cm", "Sección 1"],
       ["Diámetro zona homogénea ``D``", fmt(D_cm, 2) * " cm", "Diagonal mayor de la caja"],
       ["Campo peor caso", fmt(B_peor, 0) * " µT", "Tierra + órbita (sección 2)"],
       ["Campo objetivo por eje", fmt(B_obj, 0) * " µT", "Peor caso × factor de seguridad $(FS)"],
       ["Separación mínima ``d``", fmt(d_min, 1) * " cm", "Giroscopio + holguras + perfil (sección 3)"],
       ["Lados X / Y / Z", join(fmt.(L_ejes ./ 100, 2), " / ") * " m", "``L = d/0{,}5445``, +``\\Delta`` por par"],
       ["Separaciones X / Y / Z", join(fmt.(d_ejes, 2), " / ") * " cm", "``d = 0{,}5445\\,L``"],
       ["``NI`` X / Y / Z", join(fmt.(NI_ejes, 1), " / ") * " A·vuelta", "``B_\\text{obj}/(\\mu_0 C)`` (sección 5)"],
       ["Vueltas por bobina ``N``", string(N_vueltas), "Elección de diseño"],
       ["Corriente X / Y / Z", join(fmt.(I_ejes, 2), " / ") * " A", "``NI/N``"],
       ["Vueltas de AWG $(awg) que caben en el carril", string(N_carril) * (N_vueltas <= N_carril ? " (caben las $(N_vueltas))" : " (no caben las $(N_vueltas))"), "Carril en T del perfil (sección 5.4)"],
       ["Inductancia por par X / Y / Z", join(fmt.(1e3 .* Lpar_ejes, 1), " / ") * " mH", "``2(L_\\text{bobina} + M)`` (sección 5.5)"],
       ["Calor por par X / Y / Z", join(fmt.(Pcal_ejes, 1), " / ") * " W", "``I^2R`` en equilibrio, aire a $(T_amb) °C (sección 5.6)"],
       ["Temperatura del devanado X / Y / Z", join(fmt.(Tfin_ejes, 1), " / ") * " °C", "Convección + radiación (sección 5.6)"],
       ["Uniformidad en la zona (par $eje_ver)", "≤ " * fmt(dev_max_mag, 3) * " %", "Simulación 3D (sección 6)"],
       ["Estructura", "No ferromagnética", "Sección 7"]])

# ╔═╡ f4182795-692e-4e14-8214-da5be3b3a474
"Las 12 aristas de una caja centrada en el origen con semilados `a`, `b`, `c`."
function aristas_caja(a, b, c)
	v = [[sx * a, sy * b, sz * c] for sx in (-1, 1), sy in (-1, 1), sz in (-1, 1)]
	[(p, q) for p in v, q in v if p < q && count(p .!= q) == 1]
end

# ╔═╡ c494b369-6d32-448b-b220-cadb14eb5ad2
"Coordenadas (x, y, z) de un cuadrado de lado `L` en el plano perpendicular al eje `k`, desplazado `s` sobre ese eje."
function cuadrado(L, s, k)
	h = L / 2
	u = [h, -h, -h, h, h]
	w = [h, h, -h, -h, h]
	coords = [fill(float(s), 5), u, w]
	k == 1 && return (coords[1], coords[2], coords[3])
	k == 2 && return (coords[2], coords[1], coords[3])
	(coords[2], coords[3], coords[1])
end

# ╔═╡ 7cfeeaf7-c748-4ed6-9bc3-e17f1b9ebb60
"Superficie esférica semitransparente de radio `R` para las figuras 3D."
function esfera(R; nombre = "", color = "#2CA02C", opacidad = 0.25)
	θ = range(0, π, length = 30)
	φ = range(0, 2π, length = 50)
	surface(x = [R * sin(t) * cos(p) for t in θ, p in φ],
	        y = [R * sin(t) * sin(p) for t in θ, p in φ],
	        z = [R * cos(t) for t in θ, p in φ],
	        opacity = opacidad, showscale = false, name = nombre, showlegend = true,
	        colorscale = [[0, color], [1, color]], hoverinfo = "name")
end

# ╔═╡ 765c384d-f8c3-45fc-aeed-5e4b2fab1f02
"Ajustes comunes de las escenas 3D (misma escala en los tres ejes)."
escena3d(unidad) = attr(aspectmode = "data",
	xaxis_title = "x ($unidad)", yaxis_title = "y ($unidad)", zaxis_title = "z ($unidad)")

# ╔═╡ bc419d16-ef01-44c6-895d-08bb81c8ec53
begin
	"Bloque desplegable: muestra `titulo` y al hacer clic despliega `contenido`."
	struct Plegable{C}
		titulo::String
		contenido::C
	end
	function Base.show(io::IO, mime::MIME"text/html", p::Plegable)
		print(io, """<details style="margin: 0.6em 0; padding: 0.5em 0.9em; border: 1px solid #8884; border-radius: 6px;"><summary style="cursor: pointer; font-weight: 600;">""", p.titulo, "</summary><div style=\"padding-top: 0.5em;\">")
		show(io, mime, p.contenido)
		print(io, "</div></details>")
	end
end

# ╔═╡ 3c4d58b4-da0e-43af-8dce-d5783ddc5249
Plegable("Vocabulario mínimo", md"""
| Término | Significado |
|:---|:---|
| **Campo magnético** ``\vec B`` | La "influencia magnética" que hay en cada punto del espacio. Es un vector con magnitud y dirección. |
| **Tesla (T)** | La unidad del campo magnético. |
| **Bobina** | Alambre enrollado dando varias vueltas (**espiras**). Cuando circula corriente por él, crea un campo magnético. |
| **Corriente** ``I`` | La cantidad de carga eléctrica que pasa por el alambre cada segundo. Se mide en **amperios (A)**. |
| ``N`` | Número de vueltas de alambre de cada bobina. |
| **Homogéneo / uniforme** | El campo tiene el mismo tamaño y la misma dirección en toda una región. |
| ``\mu_0`` | Permeabilidad del vacío, una constante de la naturaleza, ``\mu_0 = 4\pi\times10^{-7}`` T·m/A.Mide la capacidad del espacio libre para permitir el paso de un campo magnetico. |
| **Ferromagnético** | Material que se imanta con facilidad: hierro, acero común, níquel, cobalto. |
""")

# ╔═╡ 2caca423-9769-488c-b5ab-b0c53e609f66
Plegable("¿De dónde sale el 0,5445?", md"""
Cerca del centro, el campo en el eje se puede escribir como una serie de potencias de ``z`` (serie de Taylor):

```math
B(z) \approx B(0) + \frac{B''(0)}{2}\,z^2 + \frac{B''''(0)}{24}\,z^4 + \dots
```

Los términos con potencias impares (``z``, ``z^3``) desaparecen por simetría: el par es igual visto desde la izquierda que desde la derecha. Si además **elegimos ``d`` para que ``B''(0) = 0``**, el primer término que deforma el campo es el de ``z^4``, que es diminuto cerca del centro. Por eso el campo queda tan plano.

Como el par son dos espiras idénticas a ``\pm d/2``, la condición ``B''(0)=0`` equivale a que **la curva de una sola espira tenga un punto de inflexión** (curvatura cero) justo en ``z = d/2``. La celda de abajo busca ese punto numéricamente con el método de bisección. El resultado es el 0,5445.
""")

# ╔═╡ 62c459b5-7d86-4409-af2e-12e907c147f8
Plegable("Deducción", md"""
1. **Un lado recto.** Para un tramo recto de alambre de largo ``2a``, visto desde un punto a distancia perpendicular ``\rho`` de su punto medio, la integral de Biot–Savart da

```math
B_\text{lado} = \frac{\mu_0 I}{4\pi\rho}\cdot\frac{2a}{\sqrt{\rho^2 + a^2}}
```

2. **Distancia de cada lado al punto del eje.** Cada lado está a ``a`` del centro de la espira y el punto está a ``z`` sobre el eje, así que ``\rho = \sqrt{a^2 + z^2}``, y ``\rho^2 + a^2 = z^2 + 2a^2``.

3. **Solo sobrevive la componente a lo largo del eje.** Cada lado produce un campo inclinado; las partes laterales de lados opuestos se cancelan. La parte que queda es ``B_\text{lado}\cdot a/\rho``.

4. **Cuatro lados iguales y ``N`` vueltas:**

```math
B = 4N\cdot\frac{\mu_0 I}{4\pi\rho}\cdot\frac{2a}{\sqrt{z^2+2a^2}}\cdot\frac{a}{\rho} = \frac{2\mu_0 N I\,a^2}{\pi\,(a^2+z^2)\sqrt{2a^2+z^2}}
```
""")

# ╔═╡ 94d0605d-447e-4c43-937f-6e85248c39f2
Plegable("¿Por qué cobre?", md"""
El alambre de cobre esmaltado es el más usado para bobinas. Su resistividad (``1{,}68\times10^{-8}\ \Omega\cdot\text{m}``) es cerca de un 37 % menor que la del aluminio (``2{,}65\times10^{-8}\ \Omega\cdot\text{m}``), así que con el mismo calibre la bobina tiene menos resistencia: la fuente necesita menos voltaje y el devanado se calienta menos. Además, el cobre se suelda y se empalma con facilidad (el aluminio se cubre de una capa de óxido que dificulta las uniones), aguanta mejor los dobleces al bobinar sin romperse y se consigue en todos los calibres.

Su desventaja es que es unas 3,3 veces más denso que el aluminio y más costoso, así que el devanado pesa más. En una jaula de laboratorio, que no se mueve, eso no es un problema.
""")

# ╔═╡ a82ff6eb-1647-4817-ac22-583f95f44b6c
Plegable("¿De dónde salen los coeficientes?", md"""
**Convección natural.** Para aire quieto y flujo laminar se usan las relaciones simplificadas de Holman (*Transferencia de calor*):

```math
h_\text{horizontal} = 1{,}32\left(\frac{\Delta T}{D}\right)^{1/4}, \qquad h_\text{vertical} = 1{,}42\left(\frac{\Delta T}{L}\right)^{1/4}
```

con ``\Delta T = T - T_\text{amb}`` en K, ``D`` el alto de la sección (perfil + lo que sobresalga del haz) y ``L`` el largo del lado, ambos en metros. Los lados verticales enfrían menos porque el aire caliente sube pegado a ellos a lo largo de todo el lado. Los pares X y Y están en planos verticales (2 lados verticales y 2 horizontales); el par Z queda acostado (4 lados horizontales).

Las fórmulas son para el nivel del mar. En Bogotá (≈ 2600 m, presión ``p \approx 75`` kPa) el aire es menos denso y enfría un poco menos, así que se multiplican por ``\sqrt{p/101{,}3\ \text{kPa}} \approx 0{,}86``.

**Radiación.** Sale de la ley de Stefan–Boltzmann, con las temperaturas en kelvin:

```math
h_\text{rad} = \varepsilon\,\sigma\left(T^2 + T_\text{amb}^2\right)\left(T + T_\text{amb}\right), \qquad \sigma = 5{,}67\times10^{-8}\ \text{W/(m}^2\text{K}^4\text{)}
```

La **emisividad** ``\varepsilon`` (entre 0 y 1) dice qué tan bien emite la superficie. El esmalte del alambre tiene ``\varepsilon \approx 0{,}9``, y la del perfil depende del acabado. Se usa el promedio de las dos, pesado por la parte del borde que ocupa cada una. Por la boca del carril se ve el alambre.

**Punto más caliente.** El calor que nace dentro del haz tiene que cruzarlo para llegar al aluminio. Un haz de alambre esmaltado sin resina conduce mal el calor (``k \approx 0{,}3`` W/(m·K)), y el aire que queda en los huecos del carril, peor todavía (``k \approx 0{,}026`` W/(m·K)). Para saber cuánto más caliente queda el interior, la sección del carril se divide en cuadritos de 0,1 mm. En cada uno, el calor que entra desde los vecinos más el que genera el alambre tiene que ser igual al que sale. Ese sistema de ecuaciones se resuelve por **diferencias finitas** con el método de Gauss–Seidel. El aluminio se toma a la temperatura del devanado, y la parte del haz que da al aire de afuera se toma aislada, que es el caso más desfavorable. El resultado es una resistencia térmica ``R_\text{haz}``, en K por cada W/m de bobina:

```math
\Delta T_\text{haz} = R_\text{haz}\,\frac{P_\text{bobina}}{4L}
```

**Tiempo de calentamiento.** Mientras no llega al equilibrio, la diferencia ``P - Q`` va calentando la masa de la bobina:

```math
C\,\frac{dT}{dt} = P(T) - Q(T), \qquad C = c_\text{Cu}\,m_\text{alambre} + c_\text{Al}\,m_\text{perfil}, \quad c_\text{Cu} = 385,\ c_\text{Al} = 897\ \text{J/(kg·K)}
```

Para el perfil se supone uno ranurado tipo 2020, que tiene un 42 % de aluminio (1,69 cm² en una sección de 2 × 2 cm). El perfil guarda mucho más calor que el alambre: pesa más, y además cada kilo de aluminio necesita más del doble de energía que uno de cobre para subir un grado. Por eso es él quien fija cuánto tarda en calentarse.
""")

# ╔═╡ 3aca9669-cfc2-47cb-9ace-58f4383e5ec6
Plegable("La fórmula del tramo recto", md"""
Para un tramo recto que va del punto ``A`` al punto ``B`` con corriente ``I``, el campo en un punto ``P`` es

```math
\vec B = \frac{\mu_0 I}{4\pi}\;\frac{|\vec a| + |\vec b|}{|\vec a|\,|\vec b|\left(|\vec a|\,|\vec b| + \vec a\cdot\vec b\right)}\;\left(\vec a\times\vec b\right),
\qquad \vec a = A - P,\quad \vec b = B - P
```

Es el resultado exacto de integrar la ley de Biot–Savart a lo largo del tramo. Se puede comprobar que sobre el eje reproduce la fórmula de la sección 4.
""")

# ╔═╡ 4e2651fe-7702-434f-8b4a-9ba19092bb8a
begin
	# PlutoPlotly 0.6.6 descarga de internet tres librerías de JavaScript en cada gráfica:
	# plotly.js (esm.sh), lodash (cdn.jsdelivr.net) e interactjs (esm.sh). Sin conexión esas
	# descargas fallan y la figura queda en blanco. `grafica` cambia esas descargas por la
	# copia de plotly.js que el paquete ya guardó en el disco y por una versión mínima de
	# las cinco funciones de lodash que usa. interactjs se importa pero no se usa.
	struct GraficaOffline
		p::PlutoPlotly.PlutoPlot
	end

	lodash_mini = """({
		isTypedArray: (o) => ArrayBuffer.isView(o) && !(o instanceof DataView),
		isPlainObject: (o) => o !== null && typeof o === "object" && [Object.prototype, null].includes(Object.getPrototypeOf(o)),
		mapValues: (o, f) => Object.fromEntries(Object.entries(o).map(([k, v]) => [k, f(v, k, o)])),
		update: (o, k, f) => (o[k] = f(o[k]), o),
		union: (...listas) => [...new Set(listas.flatMap((l) => l ?? []))],
	})"""

	function Base.show(io::IO, mime::MIME"text/html", g::GraficaOffline)
		v = PlutoPlotly.get_plotly_version()
		plotly_disco = PlutoPlotly.import_local_js(PlutoPlotly.get_local_plotly_contents(v), "default")
		html = sprint(show, mime, g.p; context = io)
		reemplazos = [
			"await import('https://cdn.jsdelivr.net/npm/lodash-es@4.17.21/+esm')" => lodash_mini,
			"await import('https://esm.sh/interactjs@1.10.19')" => "({})",
			"(await import(\"https://esm.sh/plotly.js-dist-min@$v\")).default" =>
				sprint(show, MIME"text/javascript"(), plotly_disco; context = io),
		]
		for (remoto, en_disco) in reemplazos
			occursin(remoto, html) || error("grafica: no se encontró `$remoto`. ¿Cambió la versión de PlutoPlotly?")
			html = replace(html, remoto => en_disco)
		end
		print(io, html)
	end

	"Igual que `plot(...)`, pero la figura se dibuja sin internet."
	grafica(args...; kwargs...) = GraficaOffline(plot(args...; kwargs...))
end

# ╔═╡ 52a05a54-ba81-4736-9b82-1d1b579ee5b4
let
	a, b, c = a_cm / 2, b_cm / 2, c_cm / 2
	# Aristas de la caja (NaN separa los trazos)
	X = Float64[]; Y = Float64[]; Z = Float64[]
	for (p, q) in aristas_caja(a, b, c)
		append!(X, [p[1], q[1], NaN]); append!(Y, [p[2], q[2], NaN]); append!(Z, [p[3], q[3], NaN])
	end
	trazos = [
		scatter3d(x = X, y = Y, z = Z, mode = "lines", name = "CubeSat",
		          line = attr(color = "#333", width = 5)),
		scatter3d(x = [-a, a], y = [-b, b], z = [-c, c], mode = "lines", name = "Diagonal D",
		          line = attr(color = "#2CA02C", width = 8)),
		esfera(D_cm / 2; nombre = "Esfera de diámetro D", color = "#2CA02C"),
	]
	grafica(trazos, Layout(scene = escena3d("cm"), height = 450,
	                    margin = attr(l = 0, r = 0, t = 0, b = 0)))
end

# ╔═╡ 94d5c856-20d8-4970-81a0-4eeda363e03b
let
	t = deg2rad(θ_deg)
	Tx, Ty = B_tierra, 0.0                          # flecha de la Tierra (horizontal)
	Ox, Oy = B_orbita * cos(t), B_orbita * sin(t)   # flecha deseada (órbita)
	flecha(x0, y0, x1, y1, color) = attr(x = x1, y = y1, ax = x0, ay = y0, xref = "x", yref = "y", axref = "x", ayref = "y",
		showarrow = true, arrowhead = 3, arrowsize = 1, arrowwidth = 4, arrowcolor = color, text = "")
	etiqueta(x, y, texto, color) = attr(x = x, y = y, text = texto, showarrow = false,
		font = attr(color = color, size = 14))
	lim = 70
	grafica(scatter(x = [0], y = [0], mode = "markers", marker = attr(color = "black", size = 6),
	             showlegend = false),
	     Layout(height = 380, title = "Tierra + bobinas = órbita",
	            xaxis = attr(range = [-lim, lim], title = "µT", zeroline = true),
	            yaxis = attr(range = [-10, lim], scaleanchor = "x", title = "µT"),
	            annotations = [
	                flecha(0, 0, Tx, Ty, "#4C78A8"),
	                flecha(0, 0, Ox, Oy, "#2CA02C"),
	                flecha(Tx, Ty, Ox, Oy, "#E45756"),
	                etiqueta(Tx / 2, -5, "Tierra", "#4C78A8"),
	                etiqueta(Ox / 2 - 6, Oy / 2 + 5, "Órbita (deseado)", "#2CA02C"),
	                etiqueta((Tx + Ox) / 2 + 8, (Ty + Oy) / 2 + 5, "Bobinas", "#E45756"),
	            ],
	            margin = attr(t = 40, b = 40)))
end

# ╔═╡ 7fd8a4c5-c94d-4817-9b45-c9cf61ac50be
let
	θs = 0:1:180
	Bb = [sqrt(B_tierra^2 + B_orbita^2 - 2B_tierra * B_orbita * cosd(θ)) for θ in θs]
	grafica([scatter(x = θs, y = Bb, mode = "lines", name = "|B bobinas|",
	              line = attr(color = "#E45756", width = 3)),
	      scatter(x = [θ_deg], y = [B_bobinas_θ], mode = "markers", name = "ángulo elegido",
	              marker = attr(size = 12, color = "#E45756"))],
	     Layout(height = 300, title = "Campo de las bobinas según θ",
	            xaxis_title = "Ángulo θ entre Tierra y órbita (°)",
	            yaxis_title = "µT", showlegend = false, margin = attr(t = 40, b = 50)))
end

# ╔═╡ 8274de01-5d8f-462b-838c-6f90c5b7be37
let
	L, NI = 1.0, 73.7          # bobina de 1 m y el NI de diseño del eje X
	d = r_demo * L
	zs = range(-0.8, 0.8, length = 801)
	B1 = [NI * B_espira(z - d / 2, L) * 1e6 for z in zs]   # µT
	B2 = [NI * B_espira(z + d / 2, L) * 1e6 for z in zs]
	Bt = B1 .+ B2
	zona = 100 * semiancho_uniforme(r_demo)
	grafica([scatter(x = 100zs, y = B1, name = "Bobina 1 sola", line = attr(dash = "dash", color = "#4C78A8")),
	      scatter(x = 100zs, y = B2, name = "Bobina 2 sola", line = attr(dash = "dash", color = "#F58518")),
	      scatter(x = 100zs, y = Bt, name = "Suma (lo que siente el satélite)",
	              line = attr(width = 4, color = "#E45756"))],
	     Layout(height = 420, xaxis_title = "Posición a lo largo del eje, z (cm)",
	            yaxis_title = "Campo B (µT)", legend = attr(x = 0.01, y = 0.99),
	            title = "Par de L = 1 m, NI = 73,7 A·vuelta",
	            shapes = [rect(-zona, zona, 0, 1; yref = "paper", fillcolor = "#2CA02C",
	                           opacity = 0.15, line_width = 0),
	                      vline(-50d; line_dash = "dot", line_color = "#666"),
	                      vline(50d; line_dash = "dot", line_color = "#666")],
	            margin = attr(t = 50, b = 50)))
end

# ╔═╡ b6884674-4bab-4ed5-b250-ecaabc8afa35
let
	trazos = []
	colores = ["#E45756", "#2CA02C", "#4C78A8"]
	for (k, eje) in enumerate(["X", "Y", "Z"])
		"Bobinas $eje" in mostrar || continue
		L, d = L_ejes[k], d_ejes[k]
		for s in (-1, 1)
			P = cuadrado(L, s * d / 2, k)
			push!(trazos, scatter3d(x = P[1], y = P[2], z = P[3], mode = "lines",
				name = "Par $eje (L = $(fmt(L / 100, 2)) m)", legendgroup = eje, showlegend = s == 1,
				line = attr(color = colores[k], width = 8)))
		end
		e = zeros(3); e[k] = 0.3L
		push!(trazos, scatter3d(x = [0, e[1]], y = [0, e[2]], z = [0, e[3]], mode = "lines",
			legendgroup = eje, showlegend = false, line = attr(color = colores[k], width = 4, dash = "dash")))
	end
	if "Giroscopio" in mostrar
		push!(trazos, esfera(g_cm / 2; nombre = "Giroscopio ($(fmt(g_cm, 0)) cm)", color = "#9467BD", opacidad = 0.12))
	end
	if "Zona homogénea" in mostrar
		push!(trazos, esfera(D_cm / 2; nombre = "Zona homogénea (D = $(fmt(D_cm, 1)) cm)", color = "#2CA02C", opacidad = 0.3))
	end
	if "CubeSat" in mostrar
		X = Float64[]; Y = Float64[]; Z = Float64[]
		for (p, q) in aristas_caja(a_cm / 2, b_cm / 2, c_cm / 2)
			append!(X, [p[1], q[1], NaN]); append!(Y, [p[2], q[2], NaN]); append!(Z, [p[3], q[3], NaN])
		end
		push!(trazos, scatter3d(x = X, y = Y, z = Z, mode = "lines", name = "CubeSat",
			line = attr(color = "#333", width = 4)))
	end
	if isempty(trazos)
		md"*(Activa al menos un elemento para ver la figura.)*"
	else
		grafica([trazos...], Layout(scene = escena3d("cm"), height = 600,
			margin = attr(l = 0, r = 0, t = 0, b = 0), legend = attr(x = 0, y = 1)))
	end
end

# ╔═╡ 87471ee3-8192-44e3-8493-6433a09bd89c
let
	mm(v) = 1e3 * v
	W = mm(ancho_perfil); c = W / 2
	bb, bc, p, l = mm(b_boca) / 2, mm(b_camara) / 2, mm(prof_ranura), mm(e_labio)
	# Contorno del carril en T de la cara de arriba (y = 0 es la cara exterior)
	Tx = [-bb, -bb, -bc, -bc, bc, bc, bb, bb, -bb]
	Ty = [0, -l, -l, -p, -p, -l, -l, 0, 0]
	# Las otras tres caras: el mismo carril girado alrededor del centro del perfil (0, -c),
	# si no se encima con el de arriba
	girar(x, y, g) = (cosd(g) * x - sind(g) * (y + c), sind(g) * x + cosd(g) * (y + c) - c)
	X = Float64[]; Y = Float64[]
	for g in (p + bc < c ? (0, 90, 180, 270) : (0,))
		P = girar.(Tx, Ty, g)
		append!(X, first.(P), NaN); append!(Y, last.(P), NaN)
	end
	# Espiras como círculos; las que no quedan enteras dentro del carril, en rojo
	φ = range(0, 2π, length = 25)
	r = mm(d_esmalte) / 2
	function aros(dentro)
		ax = Float64[]; ay = Float64[]
		for (x, y) in centros(capas_carril, d_esmalte)
			(y + d_esmalte / 2 <= prof_ranura + 1e-9) == dentro || continue
			append!(ax, mm(x) .+ r .* cos.(φ), NaN); append!(ay, mm(y - prof_ranura) .+ r .* sin.(φ), NaN)
		end
		ax, ay
	end
	n_dentro = min(N_vueltas, N_carril)
	trazos = [
		scatter(x = [-c, c, c, -c, -c], y = [0, 0, -W, -W, 0], mode = "lines", fill = "toself", fillcolor = "#C9CED6",
		        line = attr(color = "#8C939D"), name = "Perfil de aluminio (esquema)", hoverinfo = "skip"),
		scatter(x = X, y = Y, mode = "lines", fill = "toself", fillcolor = "white", line = attr(color = "#8C939D"),
		        showlegend = false, hoverinfo = "skip"),
		scatter(x = aros(true)[1], y = aros(true)[2], mode = "lines", fill = "toself", fillcolor = "#C87533",
		        line = attr(color = "#7A4A1E", width = 1), name = "Espira dentro del carril ($n_dentro)", hoverinfo = "skip"),
	]
	if N_vueltas > N_carril
		push!(trazos, scatter(x = aros(false)[1], y = aros(false)[2], mode = "lines", fill = "toself", fillcolor = "#E45756",
		                      line = attr(color = "#8E2A2A", width = 1), name = "Espira que no cabe ($(N_vueltas - N_carril))", hoverinfo = "skip"))
	end
	grafica(trazos, Layout(height = 480, title = "$(N_vueltas) vueltas de AWG $(awg) en el carril",
		xaxis = attr(title = "mm", range = [-c - 2, c + 2], zeroline = false, constrain = "domain"),
		yaxis = attr(title = "mm", scaleanchor = "x", range = [-W - 1, max(2, mm(sobresale) + 1.5)], zeroline = false),
		legend = attr(orientation = "h", x = 0, y = -0.15), margin = attr(t = 50, b = 40)))
end

# ╔═╡ ef9a1c2a-d36f-4ba0-ace7-9815df053870
let
	colores = ["#E45756", "#2CA02C", "#4C78A8"]
	n = min(length(t_min), 2maximum(t95_ejes) + 1)   # se muestra hasta el doble del tiempo al 95 %
	trazos = [scatter(x = t_min[1:n], y = curvas_T[k][1:n], mode = "lines", name = "Par $eje",
	                  line = attr(color = colores[k], width = 3))
	          for (k, eje) in enumerate(["X", "Y", "Z"])]
	grafica(trazos, Layout(height = 360, title = "Calentamiento de cada bobina",
	                       xaxis_title = "Tiempo (min)", yaxis_title = "Temperatura del devanado (°C)",
	                       shapes = [hline(T_amb; line_dash = "dot", line_color = "#666")],
	                       legend = attr(x = 0.99, y = 0.05, xanchor = "right"),
	                       margin = attr(t = 50, b = 50)))
end

# ╔═╡ 6d31338b-4809-4df3-b4e4-8ced1bde09d2
let
	nz, nx = 141, 121
	zs = range(-0.6L_ver, 0.6L_ver, length = nz)     # a lo largo del eje del par
	xs = range(-0.47L_ver, 0.47L_ver, length = nx)   # perpendicular al eje
	dev = [clamp(desviacion([x, 0.0, z])[1], -2.5, 2.5) for x in xs, z in zs]
	filas = [collect(dev[i, :]) for i in 1:nx]
	R = D_cm / 2
	# transpose = false: `filas` ya viene como lista de filas (PlotlyBase transpone por defecto)
	grafica([heatmap(x = 100zs, y = 100xs, z = filas, transpose = false, zmin = -2, zmax = 2, colorscale = [[0, "#2166AC"], [0.25, "#92C5DE"], [0.5, "#FFFFFF"], [0.75, "#F4A582"], [1, "#B2182B"]],
	              hovertemplate = "eje: %{x:.1f} cm<br>lateral: %{y:.1f} cm<br>desviación: %{z:.3f} %<extra></extra>",
	              colorbar = attr(title = "Desviación", tickvals = [-2, -1, 0, 1, 2],
	                              ticktext = ["−2 % (más débil)", "−1 %", "0", "+1 %", "+2 % (más fuerte)"])),
	      contour(x = 100zs, y = 100xs, z = filas, transpose = false, showscale = false, hoverinfo = "skip",
	              contours_start = -1, contours_end = 1, contours_size = 2, contours_coloring = "none",
	              line = attr(color = "black", width = 2))],
	     Layout(height = 520, title = "Par $eje_ver: L = $(fmt(L_ver, 2)) m, d = $(fmt(100d_ver, 2)) cm",
	            xaxis = attr(title = "Posición a lo largo del eje del par (cm)"),
	            yaxis = attr(title = "Posición perpendicular al eje (cm)", scaleanchor = "x"),
	            shapes = [circle(-R, R, -R, R; line_color = "#2CA02C", line_width = 4),
	                      line(-50d_ver, -50d_ver, -50L_ver, 50L_ver; line_color = "#E45756", line_width = 6),
	                      line(50d_ver, 50d_ver, -50L_ver, 50L_ver; line_color = "#E45756", line_width = 6)],
	            margin = attr(t = 50, b = 50)))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
PlutoPlotly = "~0.6.6"
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.1"
manifest_format = "2.1"
project_hash = "f37ddfae3e86173dc3f8ab4bc2068a0b0e64814b"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "e71ee7b4aa06b045259a7d6101e1cb45ad140bce"
registries = "General"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
registries = "General"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "61761f58648aa7217445f24f841839b78c712232"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.3"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
registries = "General"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Reexport"]
git-tree-sha1 = "291665b547f137df070e4dd83e432b5fee8cc4a0"
registries = "General"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.2"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.5.5+2"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
registries = "General"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
registries = "General"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
registries = "General"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
registries = "General"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
registries = "General"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
registries = "General"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "4657a834b01ce00f2b94546095c4350617b9af2c"
registries = "General"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.10.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f88f3ccef05a6a72a0cf0ed417c8fd68530f4ab2"
registries = "General"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "1.0.0"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.18.0+1"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.1+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "OpenSSL_jll", "Zlib_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.104+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.13.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
registries = "General"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2026.8.13"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "94ba93778373a53bfd5a0caaf7d809c445292ff4"
registries = "General"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.46.0+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
registries = "General"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools"]
git-tree-sha1 = "663e8b48b789916221e0765393b289ca6c88f24e"
registries = "General"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "3.0.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "Zstd_jll", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.13.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Colors", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "6256ab3ee24ef079b3afa310593817e069925eeb"
registries = "General"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.23"

    [deps.PlotlyBase.extensions]
    DataFramesExt = "DataFrames"
    DistributionsExt = "Distributions"
    IJuliaExt = "IJulia"
    JSON3Ext = "JSON3"

    [deps.PlotlyBase.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Artifacts", "ColorSchemes", "Colors", "Dates", "Downloads", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "Pkg", "PlotlyBase", "PrecompileTools", "Reexport", "ScopedValues", "Scratch", "TOML"]
git-tree-sha1 = "2b9e3d771adfe535a4fdda855f4741fdaacd3f7f"
registries = "General"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.6.6"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
registries = "General"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
registries = "General"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "5005266de4bfe50e53ff44a5cb5c540b6e47a254"
registries = "General"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.6.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.REPL]]
deps = ["Base64", "Dates", "FileWatching", "InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
registries = "General"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "67a144433c4ce877ee6d1ada69a124d6b1ecf7be"
registries = "General"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.6.2"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
registries = "General"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "e2b53ce13a53367e96601081e33d34746b571bad"
registries = "General"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.5"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "b814d5005d6a529d740ffe06f8a86396f6501138"
registries = "General"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.9.2"

    [deps.StructUtils.extensions]
    StructUtilsLazilyInitializedFieldsExt = ["LazilyInitializedFields"]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    LazilyInitializedFields = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
registries = "General"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
registries = "General"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "908fec9df6c5de98548ead82a468c95ccf6cd263"
registries = "General"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.7.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
registries = "General"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl"]
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.67.1+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.8.2+0"

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
"""

# ╔═╡ Cell order:
# ╟─3373b2f8-d052-4103-96c3-99cbc4d33f43
# ╟─80568fe3-07b0-4144-83a7-ba8700988f33
# ╠═6f80bf86-6cb4-45be-8cbe-669373fcec5b
# ╟─3c4d58b4-da0e-43af-8dce-d5783ddc5249
# ╠═41024808-f97d-41f2-9cbe-87fe50712b39
# ╟─060bd874-245b-4d00-9140-f3bc9060b16f
# ╟─93b6a733-c764-4525-a68c-1ed211664d6c
# ╟─734d8df8-56ca-4d4c-a3ef-09e737943bde
# ╟─52a05a54-ba81-4736-9b82-1d1b579ee5b4
# ╠═8ae398d3-a0cf-46b7-ae5f-807fb9d5fc6d
# ╟─4a31f1ed-1a0b-4b71-b1ad-76ed30d9d300
# ╟─79557bf3-5dc6-4a08-b56c-fedbff16e618
# ╟─94d5c856-20d8-4970-81a0-4eeda363e03b
# ╟─7fd8a4c5-c94d-4817-9b45-c9cf61ac50be
# ╟─7de83929-3370-4c88-99ab-06dfc6f9ebb6
# ╟─2654b6d3-13c5-47af-a5f8-689ca33ba2c4
# ╟─25169c00-231a-41d0-b23c-bfd9db80599e
# ╟─8274de01-5d8f-462b-838c-6f90c5b7be37
# ╟─45792165-6345-439b-ae15-d07a82184fd4
# ╟─660505fb-d06a-44e0-9d96-b5b4993510bb
# ╟─2caca423-9769-488c-b5ab-b0c53e609f66
# ╟─7fbe4978-87c0-453b-a8e5-61f67ffa9597
# ╟─bd085751-070e-4418-9499-1b19290bfcce
# ╠═df2fb722-93d5-4806-bd6a-e841537923c7
# ╟─d110feb3-2397-476c-b92f-0f2925ddacee
# ╟─e9a3333a-256c-41ef-a501-c960c758e9f3
# ╟─edfbbd1a-c924-43a1-8088-5d686074f138
# ╟─517de74b-0cf1-4f0b-8e82-643c2af149e5
# ╟─b03602d8-6c46-4061-9d46-d4cefa05a519
# ╟─0df888b1-ddc3-4436-b3b0-b20e88063bc6
# ╟─b6884674-4bab-4ed5-b250-ecaabc8afa35
# ╟─4628876c-741e-44c1-99c0-60f2785b6f5a
# ╟─62c459b5-7d86-4409-af2e-12e907c147f8
# ╟─1fda2d6f-75f7-475c-b099-c13766881b98
# ╟─fcb17b35-d8fd-4713-acb9-a5c1fe5d67dd
# ╟─29c6db60-1ed6-4c2c-b2b4-32905282cd3e
# ╟─51a44e7e-2a30-4a1f-b818-912bb15f218c
# ╟─f89b9001-46c7-4fdb-b810-117b342059b6
# ╟─0de38e59-3c0f-47e8-b716-09749ecee652
# ╟─5db17c51-8b34-4670-b558-3da87d4d4e7d
# ╟─5edc19be-f1cd-4de1-a6f8-204f033dba73
# ╠═cd93dea0-a34b-4e07-ba27-0a0a13b5ac8c
# ╟─1660aadd-b037-4477-b3cd-0b28c616d0f8
# ╟─5c827fe0-e03a-4440-9107-ef6b2d0fcab4
# ╟─e2b5b115-6061-49ec-ae58-32c77c6ccb10
# ╟─a7b47a0a-22b3-447e-885c-4b61696c72ff
# ╠═9dc5475b-fb49-478d-9477-b019e9eb0229
# ╟─1371faf2-3a2d-459b-8357-af1076cae8c3
# ╟─53e661ca-0e66-4f25-a6bb-258ce8a758b9
# ╟─94d0605d-447e-4c43-937f-6e85248c39f2
# ╟─497ac44c-20ce-41ac-9e1a-930c011639fb
# ╟─7b111d0e-e224-44e7-a660-3689a360ba09
# ╟─87471ee3-8192-44e3-8493-6433a09bd89c
# ╟─5da0dbdb-5db2-4b1c-bfa4-f841ea69f5eb
# ╠═90093a54-0f76-4ea2-8ad3-457fc5f3c93f
# ╟─32f8e23c-40ed-4653-825d-66d13c109f48
# ╟─c7bb36c0-a2ac-42a0-a1ed-a11b6dd8bfdb
# ╟─51ec3d44-bc87-4783-85f5-afa0190995b0
# ╠═478ea9ab-a3f2-472a-8d97-e4d3684d218b
# ╟─a82ff6eb-1647-4817-ac22-583f95f44b6c
# ╟─4179b755-8888-4ac1-9e90-85e98b1f90f4
# ╟─cf499b0f-54fd-4d42-9a03-a53fa848b539
# ╟─034bf3c5-e6d0-4b9c-81c6-40c1b241e749
# ╟─ef9a1c2a-d36f-4ba0-ace7-9815df053870
# ╠═14bb9727-15ca-4954-9de4-76566726a9af
# ╟─cb618cf5-b404-4d3f-af27-b74bfc1fcee8
# ╟─3aca9669-cfc2-47cb-9ace-58f4383e5ec6
# ╟─21cb32e6-3fb5-4475-8e20-575796d24dfe
# ╟─24b6b149-97cf-4d2a-8605-5fbe98a9cf29
# ╟─76baad62-fa2f-4540-9129-bded98520c8d
# ╟─75724c50-ee26-4c46-98e9-8c12ed7110fd
# ╟─025766ec-d4c4-4673-a31c-468771a2fec0
# ╟─6d31338b-4809-4df3-b4e4-8ced1bde09d2
# ╟─68f7d89d-b0ef-4a5c-88f2-e5b17a7d96b3
# ╟─0dda238e-8351-4e42-a4c9-c18c5a36a59c
# ╠═08c18c59-d31b-4aa5-b61b-d9206b1fb174
# ╟─69d3dc2e-1904-44dc-9207-a7ae8c7e995d
# ╟─845a1707-f08a-43b5-a5e8-ac71114ec5cd
# ╟─2fa60789-9ea4-471b-99aa-2b69584bb79e
# ╠═fac60950-c64c-4d88-bd27-82a60b737d8f
# ╟─127d07c1-71d9-478f-b254-6398e757c941
# ╟─ea853ecb-7073-4ffe-9b03-6b4d40d7b0ec
# ╟─96a3a5aa-183c-41c4-87d8-dc75d029e895
# ╟─70e73562-6b5f-468d-b495-0ffd677983f1
# ╟─2f65f759-0037-41f9-8eec-c4fb10dda0ed
# ╟─f4182795-692e-4e14-8214-da5be3b3a474
# ╟─c494b369-6d32-448b-b220-cadb14eb5ad2
# ╟─7cfeeaf7-c748-4ed6-9bc3-e17f1b9ebb60
# ╟─765c384d-f8c3-45fc-aeed-5e4b2fab1f02
# ╟─bc419d16-ef01-44c6-895d-08bb81c8ec53
# ╟─4e2651fe-7702-434f-8b4a-9ba19092bb8a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
