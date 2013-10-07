
def example(password):
  #the following is mostly from http://www.patrick-fuller.com/texting-from-your-server-in-python/
  import smtplib
  from email.mime.text import MIMEText

  # Message to be sent
  message = MIMEText("Hello, texting!")

  # Sending email username/password and receiving phone number
  email_username = "nickszap"
  #email_password = ""
  email_password = password
  phone_number = "4236458205"

  # Gmail to Verizon. Change here for different combinations.
  #http://en.wikipedia.org/wiki/List_of_SMS_gateways
  email_username += "@gmail.com"
  phone_number += "@txt.att.net" # for verizon "@vtext.com"

  # Format message to look like an email
  message["From"] = email_username
  message["To"] = phone_number
  message["Subject"] = "From your server!"

  # Connect and send
  s = smtplib.SMTP('smtp.gmail.com:587')
  s.starttls()
  s.login(email_username, email_password)
  s.sendmail(email_username, phone_number, message.as_string())
  s.quit()
