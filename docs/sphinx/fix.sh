grep -rl '_static' build/html | xargs sed -i 's/_static/static/g'
grep -rl '_autosummary' build/html | xargs sed -i 's/_autosummary/autosummary/g'

mv build/html/_static build/html/static
mv build/html/_autosummary build/html/autosummary

