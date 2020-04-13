# This file was generated, do not modify it. # hide
function simulate(site, n, p)
    print(site, " ")
    for i=1:n
        site = rand() < p ? rand(setdiff("ATCG", site)) : site
        print("⟶  $site ")
    end
end

simulate('A', 10, 0.4)