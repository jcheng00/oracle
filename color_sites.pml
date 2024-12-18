from pymol import cmd
cmd.select('site_0', 'resi 359-359')
cmd.color('yellow', 'site_0')
cmd.select('site_1', 'resi 384-384')
cmd.color('yellow', 'site_1')
cmd.select('site_2', 'resi 112-112')
cmd.color('blue', 'site_2')
cmd.select('site_3', 'resi 114-114')
cmd.color('blue', 'site_3')
cmd.select('site_4', 'resi 115-115')
cmd.color('blue', 'site_4')
cmd.select('site_5', 'resi 204-211')
cmd.color('blue', 'site_5')
cmd.select('site_6', 'resi 293-293')
cmd.color('blue', 'site_6')
cmd.select('site_7', 'resi 326-326')
cmd.color('blue', 'site_7')
cmd.select('site_8', 'resi 356-356')
cmd.color('blue', 'site_8')
cmd.select('site_9', 'resi 441-441')
cmd.color('blue', 'site_9')
cmd.select('site_10', 'resi 442-442')
cmd.color('blue', 'site_10')
cmd.show('sticks', 'site_*')
