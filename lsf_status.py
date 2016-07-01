from whose_running_2 import lsf_user
import subprocess

proc = subprocess.Popen(['bjobs'], stdout=subprocess.PIPE)
bjobs = proc.stdout.read()
bjobs_split = bjobs.split('\n')

jonatha = lsf_user()
queues = []
modes = []
for line in bjobs_split:
    line_split = line.split()
    if len(line_split) > 3 and line_split[3] not in queues:
        queues.append(line_split[3])
    if len(line_split) > 3 and line_split[2] not in modes:
        modes.append(line_split[2])
    jonatha.add_to_count(line_split[3], line_split[2])

print '%-9s' % '' + ''.join("%-30s" % queue for queue in queues)
# print '%-9s' % 'NAME' + ''.join(''.join("%-6s"  % mode for mode in modes) for queue in queues)