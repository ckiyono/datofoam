double reverseBytes(const char *data) // Function to reverse the byte order
{
	double result;
	char *dest = (char *)&result;
	for(unsigned int i = 0; i < sizeof(double); i++) dest[i] = data[sizeof(double)-i-1];
	return result;
}
