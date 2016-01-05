import smtplib

sender = 'jonathan.weinstein2012@gmail.com'
receivers = ['jonathan.weinstein2012@gmail.com']

message = """From: LSFManager <LSF@manager.com>
To: Me <jonathan.weinstein2012@gmail.com>
Subject: LSFManager Report

This is a test e-mail message.
"""

try:
   smtpObj = smtplib.SMTP('localhost')
   smtpObj.sendmail(sender, receivers, message)
   print("Successfully sent email")
except:
   print("Error: unable to send email")