# This file was generated, do not modify it. # hide
function simulate(site, n, p)
    print("(X(0) = $site")
    for i=1:n
        site = rand() < p ? rand(setdiff("ATCG", site)) : site
        print(", X($i) = $site")
    end
    print(")")
end

simulate('A', 5, 0.4)